"""End-to-end orchestration for AMR-Surveil.

Order of operations (each step is logged and nothing is dropped silently):

1. load + validate configuration; read tables; add canonical optional columns
2. audit the *raw* tables (quality gate; strict mode stops on ERROR findings)
3. harmonise case/whitespace label variants (recorded in ``label_harmonisation.csv``)
4. filter AMR calls by evidence tier (recorded in the exclusion log)
5. descriptive indicators, comparability flags, AST re-check against supplied breakpoints
6. data inventory, missingness, exact-ID linkage, plasmids, predictions
7. genotype-phenotype go/no-go gate, then concordance only when the gate says GO
8. write outputs, the exclusion log, a run manifest with output hashes, and the report
"""

from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from amr_surveil._version import __version__
from amr_surveil.antibiogram import cumulative_antibiogram
from amr_surveil.ast_interpretation import recheck_ast_categories
from amr_surveil.comparability import comparability_flags
from amr_surveil.config import load_config, resolve_reference_date
from amr_surveil.constants import (
    AGGREGATE_COLUMNS,
    AMR_CALL_COLUMNS,
    AST_COLUMNS,
    BREAKPOINT_COLUMNS,
    CONSUMPTION_COLUMNS,
    INCIDENCE_COLUMNS,
    PLASMID_COLUMNS,
    PREDICTION_COLUMNS,
    SAMPLE_COLUMNS,
    SCREENING_COLUMNS,
)
from amr_surveil.epidemiology import (
    aggregate_surveillance_context,
    indicators_from_frame,
    prepare_analysis_frame,
)
from amr_surveil.exceptions import QualityGateError
from amr_surveil.exclusions import ExclusionLog
from amr_surveil.genotype_phenotype import assess_readiness, genotype_phenotype_concordance
from amr_surveil.inventory import build_data_inventory, missingness_summary
from amr_surveil.io import maybe_path_hash, read_csv, sha256_file, write_csv, write_json
from amr_surveil.linkage import (
    build_matched_isolates,
    build_prediction_context,
    gene_drug_stratification,
)
from amr_surveil.manifest import unique_manifest
from amr_surveil.normalize import canonical_spellings, harmonise_labels
from amr_surveil.prediction_eval import prediction_vs_observed
from amr_surveil.reporting import build_report
from amr_surveil.surveillance import antimicrobial_consumption_context, population_incidence_context
from amr_surveil.validation import audit_tables

HARMONISED_COLUMNS = {
    "samples": ["organism", "country", "region", "source_sector", "specimen"],
    "amr_calls": ["organism"],
    "ast": ["organism", "antimicrobial", "method"],
    "predictions": ["antimicrobial"],
    "plasmids": ["organism"],
    "gene_map": ["organism", "antimicrobial"],
}

OUTPUT_FILES = [
    "audit_issues.csv",
    "exclusion_log.csv",
    "label_harmonisation.csv",
    "data_inventory.csv",
    "ast_missingness.csv",
    "resistance_by_group.csv",
    "antibiogram.csv",
    "comparability_flags.csv",
    "ast_category_recheck.csv",
    "aggregate_surveillance_context.csv",
    "incidence_context.csv",
    "consumption_context.csv",
    "matched_isolates.csv",
    "model_prediction_context.csv",
    "prediction_vs_observed.csv",
    "plasmid_calls_linked.csv",
    "gene_drug_descriptive_stratification.csv",
    "genotype_phenotype_concordance.csv",
    "genotype_phenotype_discordant.csv",
    "deduplication_log.csv",
    "readiness_assessment.json",
    "analysis_summary.json",
    "report.md",
]


def _optional(path: str | Path | None, required: list[str], table: str) -> pd.DataFrame:
    return (
        pd.DataFrame(columns=required)
        if path is None
        else read_csv(path, required=required, table=table)
    )


def _fill_columns(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    for col in columns:
        if col not in frame.columns:
            frame[col] = pd.NA
    return frame


def _filter_plasmids(
    plasmids: pd.DataFrame, samples: pd.DataFrame, excluded_tiers: set[str], log: ExclusionLog
) -> pd.DataFrame:
    """Keep plasmid rows linked to a unique manifest ID, a consistent taxon and an accepted tier."""
    if plasmids.empty:
        return plasmids.iloc[0:0].copy()
    manifest = unique_manifest(samples)
    cols = [c for c in ["sample_id", "organism"] if c in manifest.columns]
    keep_ids = set(manifest["sample_id"].astype(str))
    sid = plasmids["sample_id"].astype(str)
    for rid in sorted(set(sid[~sid.isin(keep_ids)])):
        log.add(
            "plasmids",
            "plasmids",
            rid,
            "unmatched_or_ambiguous_sample_id",
            "not linked to a unique manifest row",
        )
    linked = plasmids[sid.isin(keep_ids)].copy()
    tier = (
        linked["evidence_tier"].fillna("unspecified")
        if "evidence_tier" in linked.columns
        else pd.Series("unspecified", index=linked.index)
    )
    drop = tier.isin(excluded_tiers)
    for rid in (
        linked.loc[drop, "sample_id"].astype(str) + "|" + linked.loc[drop, "replicon"].astype(str)
    ).unique():
        log.add(
            "plasmids",
            "plasmids",
            rid,
            "excluded_evidence_tier",
            "replicon evidence tier excluded by configuration",
        )
    linked = linked[~drop]
    if "organism" in manifest.columns and not linked.empty:
        linked = linked.merge(
            manifest[cols],
            on="sample_id",
            how="inner",
            suffixes=("", "_manifest"),
            validate="many_to_one",
        )
        if "organism_manifest" in linked.columns:
            a = linked["organism"].astype("string").str.strip().str.casefold()
            b = linked["organism_manifest"].astype("string").str.strip().str.casefold()
            conflict = a.notna() & b.notna() & (a != b)
            for rid in linked.loc[conflict, "sample_id"].astype(str).unique():
                log.add(
                    "plasmids",
                    "plasmids",
                    rid,
                    "organism_mismatch",
                    "plasmid row organism conflicts with the manifest",
                )
            linked = linked[~conflict].copy()
            linked["organism"] = linked["organism"].fillna(linked["organism_manifest"])
            linked = linked.drop(columns="organism_manifest")
    return linked.reset_index(drop=True)


def _log_dropped_rows(log: ExclusionLog, table: str, raw: pd.DataFrame, kept: pd.DataFrame) -> None:
    """Record every source row that a context calculation could not use (non-strict runs)."""
    for idx in raw.index.difference(kept.index):
        log.add(
            "context_streams",
            table,
            f"row {int(idx) + 2}",
            "invalid_row_excluded_from_context",
            "row failed numeric/denominator rules; see audit_issues.csv",
        )


def _harmonise(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Fold case/whitespace variants of each label *jointly across tables*.

    Pooling matters: if the AST table spells ``escherichia coli`` and the manifest
    ``Escherichia coli``, both must end up with the same spelling or one species would split
    into two indicator groups.
    """
    records: list[dict[str, object]] = []
    by_column: dict[str, list[str]] = {}
    for name, columns in HARMONISED_COLUMNS.items():
        for col in columns:
            if col in tables[name].columns and tables[name][col].notna().any():
                by_column.setdefault(col, []).append(name)
    for col, names in by_column.items():
        canonical = canonical_spellings(pd.concat([tables[n][col] for n in names]))
        for name in names:
            tables[name][col], merged = harmonise_labels(tables[name][col], canonical)
            for item in merged:
                records.append(
                    {
                        "table": name,
                        "column": col,
                        "canonical": item["canonical"],
                        "merged_variants": "; ".join(item["variants"]),
                    }
                )
    return pd.DataFrame(records, columns=["table", "column", "canonical", "merged_variants"])


def run_pipeline(
    *,
    samples_path: str | Path,
    output_dir: str | Path,
    amr_calls_path: str | Path | None = None,
    ast_path: str | Path | None = None,
    aggregate_path: str | Path | None = None,
    predictions_path: str | Path | None = None,
    plasmids_path: str | Path | None = None,
    incidence_path: str | Path | None = None,
    consumption_path: str | Path | None = None,
    gene_drug_map_path: str | Path | None = None,
    screening_path: str | Path | None = None,
    breakpoints_path: str | Path | None = None,
    config_path: str | Path | None = None,
    config_overrides: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Run validation, descriptive analysis, exact-ID linkage and reporting.

    ``config_overrides`` (e.g. ``{"quality": {"strict": False}}``) is merged over the YAML file
    and recorded in the manifest, so exploratory runs never need a temporary config file.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path)
    if config_overrides:
        from amr_surveil.config import (  # local: avoids cycle in docs tooling
            DEFAULTS,
            validate_config,
        )

        loaded: dict[str, Any] = {k: dict(v) for k, v in config.items()}
        for section, values in config_overrides.items():
            if section not in DEFAULTS:
                validate_config({section: {}})  # raises the standard error
            loaded.setdefault(section, {}).update(values)
        config = validate_config(loaded)
    a, study = config["analysis"], config["study"]
    reference_date = resolve_reference_date(config)
    excluded_tiers = set(a["excluded_evidence_tiers"])
    log = ExclusionLog()

    # ---- 1. read
    samples = read_csv(samples_path, required=["sample_id"], table="samples")
    t: dict[str, pd.DataFrame] = {
        "samples": samples,
        "amr_calls": _optional(amr_calls_path, ["sample_id", "gene"], "amr_calls"),
        "screening": _optional(screening_path, ["sample_id"], "screening"),
        "ast": _optional(ast_path, ["sample_id", "antimicrobial", "result"], "ast"),
        "breakpoints": _optional(
            breakpoints_path, ["organism", "antimicrobial", "standard", "version"], "breakpoints"
        ),
        "aggregate": _optional(
            aggregate_path,
            [
                "year",
                "country",
                "organism",
                "antimicrobial",
                "resistant_n",
                "tested_n",
                "source_name",
                "metric_definition",
            ],
            "aggregate_surveillance",
        ),
        "predictions": _optional(predictions_path, ["sample_id", "antimicrobial"], "predictions"),
        "plasmids": _optional(plasmids_path, ["sample_id", "replicon"], "plasmids"),
        "incidence": _optional(
            incidence_path,
            [
                "year",
                "country",
                "organism",
                "antimicrobial",
                "resistant_infection_n",
                "population_n",
                "source_name",
                "case_definition",
            ],
            "incidence",
        ),
        "consumption": _optional(
            consumption_path,
            [
                "year",
                "country",
                "care_sector",
                "antimicrobial",
                "ddd_total",
                "population_n",
                "period_days",
                "source_name",
            ],
            "consumption",
        ),
        "gene_map": _optional(
            gene_drug_map_path,
            ["organism", "gene", "antimicrobial", "mapping_source"],
            "gene_drug_map",
        ),
    }
    contracts = {
        "samples": SAMPLE_COLUMNS,
        "amr_calls": AMR_CALL_COLUMNS,
        "screening": SCREENING_COLUMNS,
        "ast": AST_COLUMNS,
        "breakpoints": BREAKPOINT_COLUMNS,
        "aggregate": AGGREGATE_COLUMNS,
        "predictions": PREDICTION_COLUMNS,
        "plasmids": PLASMID_COLUMNS,
        "incidence": INCIDENCE_COLUMNS,
        "consumption": CONSUMPTION_COLUMNS,
        "gene_map": ["organism", "gene", "antimicrobial", "mapping_source", "match_type"],
    }
    for name, cols in contracts.items():
        _fill_columns(t[name], cols)
    for name in ("amr_calls", "plasmids"):  # tier labels are compared case-insensitively
        if "evidence_tier" in t[name].columns:
            t[name]["evidence_tier"] = (
                t[name]["evidence_tier"].astype("string").str.strip().str.lower()
            )
    ast = t["ast"]
    ast["result"] = ast["result"].astype("string").str.strip().str.upper()

    # ---- 2. audit raw tables, then gate
    issues = audit_tables(
        t["samples"],
        t["amr_calls"],
        t["ast"],
        t["aggregate"],
        t["predictions"],
        t["plasmids"],
        t["incidence"],
        t["consumption"],
        screening=t["screening"],
        reference_date=reference_date,
        review_identity_pct=a["review_identity_pct"],
        review_coverage_pct=a["review_coverage_pct"],
        excluded_evidence_tiers=excluded_tiers,
    )
    write_csv(issues, out / "audit_issues.csv")
    n_errors = int((issues["severity"] == "ERROR").sum()) if not issues.empty else 0
    if n_errors and config["quality"]["strict"]:
        codes = ", ".join(sorted(issues.loc[issues["severity"] == "ERROR", "code"].unique()))
        raise QualityGateError(
            f"Strict quality gate stopped the run with {n_errors} error finding(s): {codes}. "
            f"Review {out / 'audit_issues.csv'}; exploratory mode is available with quality.strict: false."
        )

    # ---- 3. harmonise labels
    harmonisation = _harmonise(t)
    write_csv(harmonisation, out / "label_harmonisation.csv")

    # ---- 3b. a record without a sample_id can never be joined: log it and remove it
    for name in ("samples", "amr_calls", "screening", "ast", "predictions", "plasmids"):
        frame_ = t[name]
        no_id = frame_["sample_id"].isna()
        if no_id.any():
            log.add_many(
                "sample_id_check",
                name,
                [f"row {int(i) + 2}" for i in frame_.index[no_id]],
                "missing_sample_id",
                "record has no sample_id and cannot be linked",
            )
            t[name] = frame_[~no_id].reset_index(drop=True)

    # ---- 4. evidence-tier filter for AMR calls
    calls_all = t["amr_calls"].dropna(subset=["sample_id", "gene"]).copy()
    calls_all["evidence_tier"] = calls_all["evidence_tier"].fillna("unspecified")
    weak_tier = calls_all["evidence_tier"].isin(excluded_tiers)
    # AMRFinderPlus PARTIAL* / INTERNAL_STOP hits are truncated or incomplete determinants: they are
    # neither a functional positive nor evidence of absence, so they are treated like weak tiers.
    non_functional = (
        calls_all["hit_method"]
        .astype("string")
        .str.upper()
        .str.match(r"^(PARTIAL|INTERNAL_STOP)")
        .fillna(False)
        .astype(bool)
    )
    drop = weak_tier | non_functional
    for rid in (
        calls_all.loc[weak_tier, "sample_id"].astype(str)
        + "|"
        + calls_all.loc[weak_tier, "gene"].astype(str)
    ).unique():
        log.add(
            "evidence_tier",
            "amr_calls",
            rid,
            "excluded_evidence_tier",
            "tier excluded from genotype analysis by configuration",
        )
    for rid in (
        calls_all.loc[non_functional & ~weak_tier, "sample_id"].astype(str)
        + "|"
        + calls_all.loc[non_functional & ~weak_tier, "gene"].astype(str)
    ).unique():
        log.add(
            "evidence_tier",
            "amr_calls",
            rid,
            "non_functional_hit_method",
            "partial or internal-stop hit: neither a positive nor a negative genotype",
        )
    calls_used = calls_all[~drop].copy()
    uncertain_calls = calls_all[drop].copy()

    # ---- 5. indicators and comparability
    frame, dedup_log = prepare_analysis_frame(
        t["samples"], t["ast"], deduplication_rule=a["deduplication_rule"], exclusions=log
    )
    indicators = indicators_from_frame(
        frame,
        confidence_level=a["confidence_level"],
        low_denominator_threshold=a["small_denominator_warning_threshold"],
    )
    flags = comparability_flags(
        frame,
        small_n=a["small_denominator_warning_threshold"],
        volume_jump_ratio=a["testing_volume_jump_ratio"],
    )
    ast_for_recheck = t["ast"].copy()
    unique_rows = unique_manifest(t["samples"])
    manifest_org = pd.Series(
        unique_rows["organism"].to_numpy(), index=unique_rows["sample_id"].astype(str)
    )
    ast_for_recheck["organism"] = ast_for_recheck["organism"].fillna(
        ast_for_recheck["sample_id"].astype(str).map(manifest_org)
    )
    recheck = recheck_ast_categories(ast_for_recheck, t["breakpoints"])

    # Genotype-phenotype analyses honour the deduplication rule, otherwise repeat isolates from
    # one patient would inflate the 2x2 counts and narrow the intervals.
    removed_ids = (
        set(dedup_log.loc[dedup_log["decision"] == "removed", "sample_id"].astype(str))
        if not dedup_log.empty
        else set()
    )
    ast_for_gp = t["ast"][~t["ast"]["sample_id"].astype(str).isin(removed_ids)].copy()

    # ---- 6. context streams and linkage
    inventory = build_data_inventory(
        t["samples"], calls_all, t["screening"], t["ast"], t["plasmids"], t["predictions"]
    )
    missingness = missingness_summary(t["samples"], t["ast"])
    aggregate_context = aggregate_surveillance_context(t["aggregate"])
    incidence_context = population_incidence_context(t["incidence"])
    consumption_context = antimicrobial_consumption_context(t["consumption"])
    _log_dropped_rows(log, "aggregate_surveillance", t["aggregate"], aggregate_context)
    _log_dropped_rows(log, "incidence", t["incidence"], incidence_context)
    _log_dropped_rows(log, "consumption", t["consumption"], consumption_context)
    matched = build_matched_isolates(t["samples"], calls_used, t["ast"])
    prediction_context = build_prediction_context(t["samples"], t["predictions"])
    pred_eval = prediction_vs_observed(
        t["samples"],
        t["predictions"],
        t["ast"],
        independence_attested=study["predictions_independence_attested"],
        confidence_level=a["confidence_level"],
        exclusions=log,
    )
    plasmids = _filter_plasmids(t["plasmids"], t["samples"], excluded_tiers, log)
    stratification = gene_drug_stratification(
        calls_used, ast_for_gp, t["gene_map"], samples=t["samples"]
    )

    # ---- 7. genotype-phenotype gate
    readiness = assess_readiness(
        t["samples"], calls_used, t["ast"], t["screening"], t["gene_map"], study
    )
    if readiness["decision"] == "GO":
        concordance, discordant = genotype_phenotype_concordance(
            t["samples"],
            calls_used,
            ast_for_gp,
            t["screening"],
            t["gene_map"],
            uncertain_calls=uncertain_calls,
            min_isolates=a["concordance_min_isolates"],
            confidence_level=a["confidence_level"],
            exclusions=log,
        )
    else:
        concordance = discordant = pd.DataFrame()
        from amr_surveil.genotype_phenotype import CONCORDANCE_COLUMNS, DISCORDANT_COLUMNS

        concordance = pd.DataFrame(columns=CONCORDANCE_COLUMNS)
        discordant = pd.DataFrame(columns=DISCORDANT_COLUMNS)

    exclusion_frame = log.to_frame()
    tables_out = {
        "exclusion_log.csv": exclusion_frame,
        "data_inventory.csv": inventory,
        "ast_missingness.csv": missingness,
        "resistance_by_group.csv": indicators,
        "antibiogram.csv": cumulative_antibiogram(indicators),
        "comparability_flags.csv": flags,
        "ast_category_recheck.csv": recheck,
        "aggregate_surveillance_context.csv": aggregate_context,
        "incidence_context.csv": incidence_context,
        "consumption_context.csv": consumption_context,
        "matched_isolates.csv": matched,
        "model_prediction_context.csv": prediction_context,
        "prediction_vs_observed.csv": pred_eval,
        "plasmid_calls_linked.csv": plasmids,
        "gene_drug_descriptive_stratification.csv": stratification,
        "genotype_phenotype_concordance.csv": concordance,
        "genotype_phenotype_discordant.csv": discordant,
        "deduplication_log.csv": dedup_log,
    }
    for name, df in tables_out.items():
        write_csv(df, out / name)
    write_json(readiness, out / "readiness_assessment.json")

    recheck_counts = recheck["status"].value_counts().to_dict() if not recheck.empty else {}
    summary: dict[str, Any] = {
        "software": "AMR-Surveil",
        "software_version": __version__,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "reference_date": reference_date.isoformat(),
        "deduplication_rule": a["deduplication_rule"],
        "confidence_level": a["confidence_level"],
        "small_denominator_warning_threshold": a["small_denominator_warning_threshold"],
        "excluded_evidence_tiers": sorted(excluded_tiers),
        "strict_quality_gate": config["quality"]["strict"],
        "study": study,
        "n_samples": int(len(t["samples"])),
        "n_amr_calls": int(len(calls_all)),
        "n_amr_calls_used": int(len(calls_used)),
        "n_screened_samples": int(t["screening"]["sample_id"].nunique())
        if not t["screening"].empty
        else 0,
        "n_ast_records": int(len(t["ast"])),
        "n_ast_indicator_groups": int(len(indicators)),
        "n_low_denominator_groups": int(indicators["small_denominator_warning"].sum())
        if not indicators.empty
        else 0,
        "n_comparability_flags": int(len(flags)),
        "n_aggregate_rows": int(len(aggregate_context)),
        "n_incidence_rows": int(len(incidence_context)),
        "n_consumption_rows": int(len(consumption_context)),
        "n_predictions": int(len(t["predictions"])),
        "n_plasmid_calls": int(len(plasmids)),
        "n_matched_isolates": int(len(matched)),
        "n_gene_drug_groups": int(len(stratification)),
        "n_concordance_strata": int(len(concordance)),
        "n_discordant_records": int(len(discordant)),
        "ast_recheck_status_counts": {str(k): int(v) for k, v in recheck_counts.items()},
        "n_excluded_records": int(len(exclusion_frame)),
        "genotype_phenotype_decision": readiness["decision"],
        "n_qa_errors": n_errors,
        "n_qa_warnings": int((issues["severity"] == "WARNING").sum()) if not issues.empty else 0,
        "n_qa_info": int((issues["severity"] == "INFO").sum()) if not issues.empty else 0,
        "n_label_harmonisations": int(len(harmonisation)),
        "analysis_boundaries": [
            "Aggregate resistance, population incidence and antimicrobial consumption streams are never joined to isolate-level genomes by country, year or organism.",
            "DDD consumption is a standardised population indicator, not patient-level ingestion or individual exposure.",
            "Genotype-phenotype joins use exact sample_id only; no fuzzy matching is used.",
            "A negative genotype is defined only for isolates declared as screened; a missing call row is never read as a negative.",
            "Genotype-phenotype concordance runs only when the go/no-go gate passes and is descriptive, not diagnostic accuracy.",
            "Model predictions remain predictions: never converted to observed AST, and accuracy-type metrics appear only when independence from training data is attested.",
            "Genomic relatedness is not evidence of transmission; no transmission, diagnostic, treatment or clinical susceptibility claim is produced.",
        ],
    }
    write_json(summary, out / "analysis_summary.json")
    build_report(
        summary=summary,
        indicators=indicators,
        flags=flags,
        issues=issues,
        readiness=readiness,
        concordance=concordance,
        discordant=discordant,
        recheck=recheck,
        missingness=missingness,
        exclusions=log.summary(),
        harmonisation=harmonisation,
        pred_eval=pred_eval,
        aggregates=aggregate_context,
        incidence=incidence_context,
        consumption=consumption_context,
        output_path=out / "report.md",
    )

    # ---- 8. manifest last, so it can hash every other output
    outputs = {}
    for name in OUTPUT_FILES:
        p = out / name
        if p.is_file():
            outputs[name] = sha256_file(p)
    inputs = {
        "samples": samples_path,
        "amr_calls": amr_calls_path,
        "screening": screening_path,
        "ast": ast_path,
        "breakpoints": breakpoints_path,
        "aggregate_surveillance": aggregate_path,
        "predictions": predictions_path,
        "plasmids": plasmids_path,
        "incidence": incidence_path,
        "consumption": consumption_path,
        "gene_drug_map": gene_drug_map_path,
        "config": config_path,
    }
    write_json(
        {
            "manifest_schema_version": 1,
            "software": "AMR-Surveil",
            "software_version": __version__,
            "created_at_utc": summary["created_at_utc"],
            "environment": {
                "python": sys.version.split()[0],
                "pandas": pd.__version__,
                "platform": platform.platform(),
            },
            "effective_configuration": config,
            "inputs": {k: maybe_path_hash(v) for k, v in inputs.items()},
            "outputs_sha256": outputs,
            "limitations": summary["analysis_boundaries"],
        },
        out / "run_manifest.json",
    )
    return summary
