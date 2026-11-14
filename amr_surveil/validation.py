"""Surveillance-specific data-quality checks. Findings are never silently repaired.

Severity convention
-------------------
``ERROR``    ambiguous or impossible data that would make a join or count unreliable; stops a
             strict run.
``WARNING``  data that can be analysed but needs scientific review before interpretation.
``INFO``     provenance gaps or facts worth recording; never block a run.

Row indices in ``examples`` are 0-based data-row positions in the input file (the header is
not counted).
"""

from __future__ import annotations

import datetime as dt
import json
from collections.abc import Iterable

import numpy as np
import pandas as pd

from amr_surveil.constants import (
    DEFAULT_EXCLUDED_EVIDENCE_TIERS,
    EVIDENCE_TIERS,
    ISSUE_COLUMNS,
    VALID_AST_CATEGORIES,
    VALID_CARE_SECTORS,
    VALID_CONSUMPTION_SOURCES,
    VALID_SECTORS,
)
from amr_surveil.normalize import harmonise_labels, label_key, parse_iso_dates, parse_mic_series

_SEVERITY_ORDER = {"ERROR": 0, "WARNING": 1, "INFO": 2}


def issue_row(
    table: str, severity: str, code: str, message: str, values: Iterable[object] | None = None
) -> dict[str, object]:
    """Build one audit finding."""
    items = list(values) if values is not None else []
    return {
        "table": table,
        "severity": severity,
        "code": code,
        "message": message,
        "row_count": len(items),
        "examples": json.dumps([str(x) for x in items[:5]], ensure_ascii=False),
    }


class _Findings:
    def __init__(self) -> None:
        self.rows: list[dict[str, object]] = []

    def add(
        self,
        table: str,
        severity: str,
        code: str,
        message: str,
        values: Iterable[object] | None = None,
    ) -> None:
        self.rows.append(issue_row(table, severity, code, message, values))

    def frame(self) -> pd.DataFrame:
        frame = pd.DataFrame(self.rows, columns=ISSUE_COLUMNS)
        if frame.empty:
            return frame
        frame["_rank"] = frame["severity"].map(_SEVERITY_ORDER)
        return (
            frame.sort_values("_rank", kind="stable").drop(columns="_rank").reset_index(drop=True)
        )


def _idx(frame: pd.DataFrame, mask: pd.Series) -> list[int]:
    return frame.index[mask.fillna(False).astype(bool)].tolist()


def _label_variants(findings: _Findings, table: str, frame: pd.DataFrame, column: str) -> None:
    if column not in frame.columns or frame[column].dropna().empty:
        return
    _, merged = harmonise_labels(frame[column])
    if merged:
        examples = [" | ".join(m["variants"]) for m in merged]  # type: ignore[arg-type]
        findings.add(
            table,
            "WARNING",
            f"inconsistent_{column}_labels",
            f"Labels in '{column}' differ only by case/whitespace; they are merged onto the most frequent spelling for "
            "analysis. Synonyms (e.g. abbreviations) are not merged: curate a controlled vocabulary upstream.",
            examples,
        )


def _duplicate_key_check(
    findings: _Findings, table: str, frame: pd.DataFrame, key: list[str], code: str, what: str
) -> None:
    if frame.empty:
        return
    # An absent optional key column behaves like an all-missing one, so a caller who omits it
    # still gets the duplicate finding.
    work = frame.reindex(columns=key)
    folded = work.apply(lambda col: col.map(label_key)).astype("object").where(work.notna(), "<NA>")
    dup = folded.duplicated(keep=False)
    if dup.any():
        findings.add(
            table,
            "ERROR",
            code,
            f"The same {what} key appears more than once. Rows are never summed automatically because duplicates "
            "would double-count; remove the duplicate or add a field that distinguishes the rows.",
            _idx(frame, dup),
        )


def audit_tables(
    samples: pd.DataFrame,
    amr_calls: pd.DataFrame,
    ast: pd.DataFrame,
    aggregate: pd.DataFrame,
    predictions: pd.DataFrame,
    plasmids: pd.DataFrame,
    incidence: pd.DataFrame | None = None,
    consumption: pd.DataFrame | None = None,
    *,
    screening: pd.DataFrame | None = None,
    reference_date: dt.date | None = None,
    review_identity_pct: float = 90.0,
    review_coverage_pct: float = 60.0,
    excluded_evidence_tiers: Iterable[str] = DEFAULT_EXCLUDED_EVIDENCE_TIERS,
) -> pd.DataFrame:
    """Audit isolate, genomic, AST, surveillance, consumption and prediction streams."""
    f = _Findings()
    incidence = incidence if incidence is not None else pd.DataFrame()
    consumption = consumption if consumption is not None else pd.DataFrame()
    screening = screening if screening is not None else pd.DataFrame()
    reference = reference_date or dt.datetime.now(dt.timezone.utc).date()
    excluded_tiers = set(excluded_evidence_tiers)

    _audit_samples(f, samples, reference)
    known_ids = set(samples.get("sample_id", pd.Series(dtype=str)).dropna().astype(str))
    _audit_cross_table_ids(f, samples, known_ids, amr_calls, ast, predictions, plasmids, screening)
    _audit_amr_calls(
        f, amr_calls, screening, review_identity_pct, review_coverage_pct, excluded_tiers
    )
    _audit_screening(f, screening)
    _audit_ast(f, ast, samples)
    _audit_aggregate(f, aggregate)
    _audit_predictions(f, predictions)
    _audit_incidence(f, incidence)
    _audit_consumption(f, consumption)
    return f.frame()


# ------------------------------------------------------------------------ samples
def _audit_samples(f: _Findings, samples: pd.DataFrame, reference: dt.date) -> None:
    if samples.empty:
        f.add("samples", "ERROR", "empty_samples", "The isolate/sample manifest contains no rows.")
        return
    if "sample_id" in samples.columns:
        blank = _idx(samples, samples["sample_id"].isna())
        if blank:
            f.add(
                "samples",
                "ERROR",
                "missing_sample_id",
                "Sample identifiers must be populated.",
                blank,
            )
        dup_mask = samples["sample_id"].notna() & samples["sample_id"].duplicated(keep=False)
        if dup_mask.any():
            f.add(
                "samples",
                "ERROR",
                "duplicate_sample_id",
                "sample_id must be unique in samples.csv; joins would be ambiguous.",
                samples.loc[dup_mask, "sample_id"].tolist(),
            )
    for col in ("organism", "collection_date", "source_sector"):
        if col in samples.columns:
            missing = _idx(samples, samples[col].isna())
            if missing:
                f.add(
                    "samples",
                    "WARNING",
                    f"missing_{col}",
                    f"Missing {col} reduces stratification or interpretability.",
                    missing,
                )
    for col in ("organism", "country", "source_sector", "specimen"):
        _label_variants(f, "samples", samples, col)

    if "collection_date" in samples.columns:
        dates = parse_iso_dates(samples["collection_date"])
        if dates["non_iso"].any():
            f.add(
                "samples",
                "WARNING",
                "non_iso_collection_date",
                "Collection dates must be ISO 8601 (YYYY, YYYY-MM or YYYY-MM-DD). Other forms such as 05/06/2024 are "
                "ambiguous (day/month order) and are not parsed; the affected year is unknown.",
                samples.loc[dates["non_iso"], "collection_date"].tolist(),
            )
        if dates["impossible"].any():
            f.add(
                "samples",
                "WARNING",
                "invalid_collection_date",
                "Collection dates are in ISO form but are not real calendar dates.",
                samples.loc[dates["impossible"], "collection_date"].tolist(),
            )
        future = dates["date"].notna() & (dates["date"] > pd.Timestamp(reference))
        if future.any():
            f.add(
                "samples",
                "ERROR",
                "future_collection_date",
                f"Collection date is after the reference date {reference.isoformat()}.",
                samples.loc[future, "collection_date"].tolist(),
            )
    if "source_sector" in samples.columns:
        sectors = samples["source_sector"].dropna().astype(str).str.strip().str.lower()
        unknown = sectors[~sectors.isin(VALID_SECTORS)].unique().tolist()
        if unknown:
            f.add(
                "samples",
                "WARNING",
                "nonstandard_source_sector",
                "Unrecognised source-sector labels are retained but should be mapped to the controlled vocabulary.",
                unknown,
            )
    for col in ("sequencing_accession", "source_study", "data_license"):
        if col in samples.columns:
            missing = _idx(samples, samples[col].isna())
            if missing:
                f.add(
                    "samples",
                    "INFO",
                    f"missing_{col}",
                    f"{col} is not recorded; provenance/licence cannot be traced for these rows.",
                    missing,
                )


def _audit_cross_table_ids(
    f: _Findings,
    samples: pd.DataFrame,
    known_ids: set[str],
    amr_calls,
    ast,
    predictions,
    plasmids,
    screening,
) -> None:
    for table_name, frame in (
        ("amr_calls", amr_calls),
        ("ast", ast),
        ("predictions", predictions),
        ("plasmids", plasmids),
        ("screening", screening),
    ):
        if frame.empty or "sample_id" not in frame.columns:
            continue
        missing_ids = _idx(frame, frame["sample_id"].isna())
        if missing_ids:
            f.add(
                table_name,
                "WARNING",
                "missing_sample_id",
                "Rows without sample_id cannot be linked to the sample manifest.",
                missing_ids,
            )
        orphaned = sorted(set(frame["sample_id"].dropna().astype(str)) - known_ids)
        if orphaned:
            f.add(
                table_name,
                "WARNING",
                "unmatched_sample_id",
                "Rows reference sample IDs absent from samples.csv; they are excluded from isolate-linked outputs.",
                orphaned,
            )
        if "organism" in frame.columns and "organism" in samples.columns:
            left = frame[["sample_id", "organism"]].dropna().drop_duplicates()
            right = samples[["sample_id", "organism"]].dropna().drop_duplicates()
            compared = left.merge(right, on="sample_id", suffixes=("_input", "_manifest"))
            mismatch = (
                compared.loc[
                    compared["organism_input"].map(label_key)
                    != compared["organism_manifest"].map(label_key),
                    "sample_id",
                ]
                .astype(str)
                .tolist()
            )
            if mismatch:
                f.add(
                    table_name,
                    "WARNING",
                    "organism_mismatch",
                    "The organism label conflicts with the sample manifest; affected isolates are excluded from paired evidence or AST indicators.",
                    mismatch,
                )


# ---------------------------------------------------------------------- amr calls
def _audit_amr_calls(
    f, amr_calls, screening, review_identity, review_coverage, excluded_tiers
) -> None:
    if amr_calls.empty:
        return
    for col in ("sample_id", "gene"):
        if col in amr_calls.columns:
            missing = _idx(amr_calls, amr_calls[col].isna())
            if missing:
                f.add(
                    "amr_calls",
                    "WARNING",
                    f"missing_{col}",
                    f"AMR call rows without {col} cannot support isolate-linked interpretation.",
                    missing,
                )
    numeric: dict[str, pd.Series] = {}
    for col in ("identity_pct", "coverage_pct"):
        if col in amr_calls.columns:
            values = pd.to_numeric(amr_calls[col], errors="coerce")
            numeric[col] = values
            bad = _idx(
                amr_calls, amr_calls[col].notna() & (values.isna() | (values < 0) | (values > 100))
            )
            if bad:
                f.add(
                    "amr_calls",
                    "ERROR",
                    f"invalid_{col}",
                    f"{col} must be numeric and between 0 and 100 when supplied.",
                    bad,
                )
    low = pd.Series(False, index=amr_calls.index)
    if "identity_pct" in numeric:
        low |= numeric["identity_pct"] < review_identity
    if "coverage_pct" in numeric:
        low |= numeric["coverage_pct"] < review_coverage
    if low.any():
        f.add(
            "amr_calls",
            "WARNING",
            "low_identity_or_coverage_review",
            f"Calls below identity {review_identity:g}% or coverage {review_coverage:g}% are flagged for expert review. "
            "These are review flags, not call filters; no call is dropped for this reason.",
            _idx(amr_calls, low),
        )
    for col in ("method", "database", "database_version"):
        if col in amr_calls.columns:
            missing = _idx(amr_calls, amr_calls[col].isna())
            if missing:
                f.add(
                    "amr_calls",
                    "WARNING",
                    f"missing_call_{col}",
                    f"AMR calls without {col} cannot be reproduced; record the tool, database and version.",
                    missing,
                )
    if {"sample_id", "gene", "method", "database", "database_version"}.issubset(amr_calls.columns):
        key = ["sample_id", "gene", "method", "database", "database_version"]
        if "contig" in amr_calls.columns:
            key += ["contig", "start", "end"]
        dup = amr_calls.duplicated(key, keep=False)
        if dup.any():
            f.add(
                "amr_calls",
                "WARNING",
                "duplicate_amr_call",
                "Repeated calls with the same sample, gene, location, method and database version should be reviewed.",
                amr_calls.loc[dup, "sample_id"].tolist(),
            )
    if "evidence_tier" in amr_calls.columns:
        tiers = amr_calls["evidence_tier"].astype("string").str.strip().str.lower()
        invalid = tiers.notna() & ~tiers.isin(EVIDENCE_TIERS)
        if invalid.any():
            f.add(
                "amr_calls",
                "ERROR",
                "invalid_evidence_tier",
                f"evidence_tier must be one of: {', '.join(EVIDENCE_TIERS)}.",
                tiers[invalid].unique().tolist(),
            )
        weak = tiers.isin(excluded_tiers)
        if weak.any():
            f.add(
                "amr_calls",
                "WARNING",
                "low_confidence_evidence_tier",
                "Calls from heuristic-motif or generic homology-search methods are not validated resistance-gene detections. "
                "They are excluded from linkage and concordance by the configured evidence-tier policy.",
                _idx(amr_calls, weak),
            )
        unspecified = tiers.isna()
        if unspecified.any():
            f.add(
                "amr_calls",
                "INFO",
                "evidence_tier_unspecified",
                "No evidence tier was supplied; these calls are retained and treated as unspecified.",
                _idx(amr_calls, unspecified),
            )
    if "hit_method" in amr_calls.columns:
        hit = amr_calls["hit_method"].astype("string").str.upper()
        partial = hit.str.startswith("PARTIAL").fillna(False) | hit.str.startswith(
            "INTERNAL_STOP"
        ).fillna(False)
        if partial.any():
            f.add(
                "amr_calls",
                "WARNING",
                "partial_or_internal_stop_hit",
                "Partial or internal-stop hits may not encode a functional determinant; review gene integrity before attributing resistance.",
                _idx(amr_calls, partial),
            )
    if (
        not screening.empty
        and "sample_id" in screening.columns
        and "sample_id" in amr_calls.columns
    ):
        unscreened = sorted(
            set(amr_calls["sample_id"].dropna().astype(str))
            - set(screening["sample_id"].dropna().astype(str))
        )
        if unscreened:
            f.add(
                "amr_calls",
                "WARNING",
                "calls_for_unscreened_sample",
                "Calls exist for samples absent from the screening table; their negative status cannot be established.",
                unscreened,
            )


def _audit_screening(f: _Findings, screening: pd.DataFrame) -> None:
    if screening.empty:
        return
    for col in ("method", "database", "database_version"):
        if col in screening.columns:
            missing = _idx(screening, screening[col].isna())
            if missing:
                f.add(
                    "screening",
                    "WARNING",
                    f"missing_screening_{col}",
                    f"Screening records without {col} are not reproducible.",
                    missing,
                )
    cols = [c for c in ("method", "database", "database_version") if c in screening.columns]
    if "sample_id" in screening.columns and cols:
        per_sample = (
            screening.dropna(subset=["sample_id"])
            .groupby("sample_id")[cols]
            .apply(lambda g: len(g.drop_duplicates()))
        )
        multi = per_sample[per_sample > 1]
        if not multi.empty:
            f.add(
                "screening",
                "WARNING",
                "multiple_screening_records",
                "A sample has more than one distinct screening method/database; it is excluded from concordance because the negative status is ambiguous.",
                multi.index.astype(str).tolist(),
            )


# ---------------------------------------------------------------------------- AST
def _audit_ast(f: _Findings, ast: pd.DataFrame, samples: pd.DataFrame) -> None:
    if ast.empty:
        return
    for col in ("antimicrobial", "result"):
        if col in ast.columns:
            missing = _idx(ast, ast[col].isna())
            if missing:
                f.add(
                    "ast",
                    "WARNING",
                    f"missing_{col}",
                    f"AST rows without {col} cannot contribute to categorical indicators.",
                    missing,
                )
    _label_variants(f, "ast", ast, "antimicrobial")
    if "result" in ast.columns:
        result = ast["result"].astype("string").str.strip().str.upper()
        invalid = ast["result"].notna() & ~result.isin(VALID_AST_CATEGORIES)
        if invalid.any():
            f.add(
                "ast",
                "ERROR",
                "invalid_ast_category",
                "AST result must be S, I or R. EUCAST I means susceptible, increased exposure; other standards may define categories differently.",
                _idx(ast, invalid),
            )
    for col in ("method", "breakpoint_standard", "breakpoint_version"):
        if col in ast.columns:
            missing = _idx(ast, ast[col].isna())
            if missing:
                f.add(
                    "ast",
                    "WARNING",
                    f"missing_{col}",
                    f"Missing {col} prevents confident breakpoint-context interpretation.",
                    missing,
                )
    if "mic" in ast.columns:
        parsed = parse_mic_series(ast["mic"])
        present = ast["mic"].notna()
        bad = present & ~parsed["mic_valid"]
        if bad.any():
            f.add(
                "ast",
                "WARNING",
                "invalid_mic",
                "MIC must be a number, optionally with a censoring operator (<=, <, >, >=, =); decimal commas are not accepted. "
                "Original rows are retained for review.",
                ast.loc[bad, "mic"].tolist(),
            )
        nonpositive = present & parsed["mic_valid"] & (parsed["mic_value"] <= 0)
        if nonpositive.any():
            f.add(
                "ast",
                "WARNING",
                "nonpositive_mic",
                "MIC values must be greater than zero.",
                _idx(ast, nonpositive),
            )
        censored = (
            present
            & parsed["mic_valid"]
            & (parsed["mic_operator"] != "=")
            & ~ast["mic"].astype("string").str.fullmatch(r"\s*[0-9.eE+-]+\s*").fillna(False)
        )
        if censored.any():
            f.add(
                "ast",
                "INFO",
                "censored_mic_present",
                "Some MIC values are censored (<=, <, >, >=); they are never treated as exact values.",
                _idx(ast, censored),
            )
    if "zone_mm" in ast.columns:
        zone = pd.to_numeric(ast["zone_mm"], errors="coerce")
        bad_zone = ast["zone_mm"].notna() & (zone.isna() | (zone <= 0) | ~np.isfinite(zone))
        if bad_zone.any():
            f.add(
                "ast",
                "WARNING",
                "invalid_zone_mm",
                "Disk-diffusion zone diameters must be positive numbers (mm).",
                _idx(ast, bad_zone),
            )
    if {
        "sample_id",
        "antimicrobial",
        "method",
        "breakpoint_standard",
        "breakpoint_version",
        "result",
    }.issubset(ast.columns):
        key = ["sample_id", "method", "breakpoint_standard", "breakpoint_version"]
        work = ast.assign(
            _drug=ast["antimicrobial"].map(label_key),
            _res=ast["result"].astype("string").str.strip().str.upper(),
        )
        gkey = key + ["_drug"]
        counts = work.groupby(gkey, dropna=False)["_res"].nunique()
        conflicts = counts[counts > 1]
        if not conflicts.empty:
            f.add(
                "ast",
                "ERROR",
                "conflicting_duplicate_ast",
                "Conflicting S/I/R records share the same isolate/drug/method/breakpoint key; resolve before aggregation.",
                [tuple(i) for i in conflicts.index[:5]]
                + ([] if len(conflicts) <= 5 else ["…"] * (len(conflicts) - 5)),
            )
        exact = work[work["_res"].notna()].duplicated(gkey + ["_res"], keep=False)
        if exact.any():
            f.add(
                "ast",
                "WARNING",
                "duplicate_ast_record",
                "Repeated classifications with the same isolate/drug/method/breakpoint/category are counted once in indicators; review differing MIC/date metadata.",
                work.loc[exact[exact].index, "sample_id"].tolist(),
            )
    if "interpretation_date" in ast.columns:
        interp = parse_iso_dates(ast["interpretation_date"])
        partial = interp["date"].notna() & (interp["precision"] != "day")
        interp["non_iso"] = interp["non_iso"] | partial
        if interp["non_iso"].any():
            f.add(
                "ast",
                "WARNING",
                "non_iso_interpretation_date",
                "interpretation_date must be ISO 8601 (YYYY-MM-DD).",
                ast.loc[interp["non_iso"], "interpretation_date"].tolist(),
            )
        if "collection_date" in samples.columns and "sample_id" in samples.columns:
            collected = samples.drop_duplicates("sample_id", keep=False)[
                ["sample_id", "collection_date"]
            ].copy()
            collected["_collected"] = parse_iso_dates(collected["collection_date"])["date"]
            joined = ast.assign(_interp=interp["date"].where(~partial)).merge(
                collected[["sample_id", "_collected"]], on="sample_id", how="left"
            )
            before = (
                joined["_interp"].notna()
                & joined["_collected"].notna()
                & (joined["_interp"] < joined["_collected"])
            )
            if before.any():
                f.add(
                    "ast",
                    "WARNING",
                    "interpretation_before_collection",
                    "An AST result is dated before its isolate's collection date.",
                    joined.loc[before, "sample_id"].tolist(),
                )


# ---------------------------------------------------- aggregate / incidence / DDD
def _audit_years(
    f: _Findings, table: str, frame: pd.DataFrame, *, include_malformed: bool = True
) -> None:
    """A surveillance year must be a whole calendar year in a plausible range."""
    if "year" not in frame.columns:
        return
    years = pd.to_numeric(frame["year"], errors="coerce")
    out_of_range = years.notna() & ~years.between(1900, 2100)
    if include_malformed:
        out_of_range = out_of_range | (frame["year"].notna() & (years.isna() | (years % 1 != 0)))
    bad = _idx(frame, out_of_range)
    if bad:
        f.add(
            table,
            "ERROR",
            "implausible_year",
            "year must be a whole calendar year between 1900 and 2100.",
            bad,
        )


def _audit_aggregate(f: _Findings, aggregate: pd.DataFrame) -> None:
    if aggregate.empty:
        return
    if "year" in aggregate.columns:
        years = pd.to_numeric(aggregate["year"], errors="coerce")
        bad_years = _idx(aggregate, aggregate["year"].notna() & (years.isna() | (years % 1 != 0)))
        if bad_years:
            f.add(
                "aggregate_surveillance",
                "WARNING",
                "invalid_aggregate_year",
                "Aggregate surveillance year should be a whole-number calendar year.",
                bad_years,
            )
    _audit_years(f, "aggregate_surveillance", aggregate, include_malformed=False)
    if {"resistant_n", "tested_n"}.issubset(aggregate.columns):
        resistant = pd.to_numeric(aggregate["resistant_n"], errors="coerce")
        tested = pd.to_numeric(aggregate["tested_n"], errors="coerce")
        bad = _idx(
            aggregate,
            resistant.isna()
            | tested.isna()
            | (resistant < 0)
            | (tested < 0)
            | (resistant > tested),
        )
        if bad:
            f.add(
                "aggregate_surveillance",
                "ERROR",
                "invalid_aggregate_counts",
                "Counts must be numeric, non-negative, and resistant_n must not exceed tested_n.",
                bad,
            )
        fractional = _idx(
            aggregate,
            resistant.notna() & tested.notna() & ((resistant % 1 != 0) | (tested % 1 != 0)),
        )
        if fractional:
            f.add(
                "aggregate_surveillance",
                "ERROR",
                "non_integer_aggregate_counts",
                "Surveillance numerator and denominator counts must be whole numbers.",
                fractional,
            )
        zero = _idx(aggregate, tested.notna() & (tested == 0))
        if zero:
            f.add(
                "aggregate_surveillance",
                "WARNING",
                "zero_aggregate_denominator",
                "Rows with tested_n=0 cannot yield a resistance proportion and are omitted from the percentage output.",
                zero,
            )
    for col in ("source_name", "metric_definition"):
        if col in aggregate.columns:
            missing = _idx(aggregate, aggregate[col].isna())
            if missing:
                f.add(
                    "aggregate_surveillance",
                    "WARNING",
                    f"missing_{col}",
                    f"Aggregate results without {col} have weak provenance.",
                    missing,
                )
    _duplicate_key_check(
        f,
        "aggregate_surveillance",
        aggregate,
        [
            "year",
            "country",
            "source_sector",
            "organism",
            "antimicrobial",
            "source_name",
            "metric_definition",
        ],
        "duplicate_aggregate_key",
        "aggregate-surveillance",
    )


def _audit_predictions(f: _Findings, predictions: pd.DataFrame) -> None:
    if predictions.empty:
        return
    if {"predicted_mic", "predicted_class"}.issubset(predictions.columns):
        missing_output = _idx(
            predictions, predictions["predicted_mic"].isna() & predictions["predicted_class"].isna()
        )
        if missing_output:
            f.add(
                "predictions",
                "WARNING",
                "missing_prediction_value",
                "Prediction rows need a predicted_mic or predicted_class value.",
                missing_output,
            )
    if "prediction_score" in predictions.columns:
        scores = pd.to_numeric(predictions["prediction_score"], errors="coerce")
        bad = _idx(
            predictions,
            predictions["prediction_score"].notna() & (scores.isna() | (scores < 0) | (scores > 1)),
        )
        if bad:
            f.add(
                "predictions",
                "WARNING",
                "non_probability_prediction_score",
                "prediction_score is expected to be in [0, 1]; inspect the model's score semantics or map it to a different field.",
                bad,
            )
    for col in ("model_version", "checkpoint_id"):
        if col in predictions.columns:
            missing = _idx(predictions, predictions[col].isna())
            if missing:
                f.add(
                    "predictions",
                    "INFO",
                    f"missing_{col}",
                    f"Predictions without {col} cannot be tied to a specific model artefact.",
                    missing,
                )
    _label_variants(f, "predictions", predictions, "antimicrobial")


def _audit_incidence(f: _Findings, incidence: pd.DataFrame) -> None:
    if incidence.empty:
        return
    for col in ("resistant_infection_n", "population_n"):
        if col not in incidence.columns:
            continue
        values = pd.to_numeric(incidence[col], errors="coerce")
        missing = _idx(incidence, values.isna())
        if missing:
            f.add(
                "incidence",
                "ERROR",
                f"missing_{col}",
                f"{col} is required for population-based incidence.",
                missing,
            )
        negative = _idx(incidence, values.notna() & (values < 0))
        if negative:
            f.add("incidence", "ERROR", f"negative_{col}", f"{col} cannot be negative.", negative)
        infinite = _idx(incidence, values.notna() & ~np.isfinite(values))
        if infinite:
            f.add("incidence", "ERROR", f"non_finite_{col}", f"{col} must be finite.", infinite)
        fractional = _idx(incidence, values.notna() & (values % 1 != 0))
        if fractional:
            f.add(
                "incidence",
                "ERROR",
                f"non_integer_{col}",
                f"{col} must be a whole-number count.",
                fractional,
            )
    if "population_multiplier" in incidence.columns:
        raw = incidence["population_multiplier"]
        parsed = pd.to_numeric(raw, errors="coerce")
        invalid_text = _idx(incidence, raw.notna() & parsed.isna())
        if invalid_text:
            f.add(
                "incidence",
                "ERROR",
                "invalid_population_multiplier",
                "population_multiplier must be numeric when supplied.",
                invalid_text,
            )
        non_positive = _idx(incidence, parsed.notna() & ((parsed <= 0) | ~np.isfinite(parsed)))
        if non_positive:
            f.add(
                "incidence",
                "ERROR",
                "invalid_population_multiplier",
                "population_multiplier must be positive.",
                non_positive,
            )
    _audit_years(f, "incidence", incidence)
    if {"resistant_infection_n", "population_n"}.issubset(incidence.columns):
        cases = pd.to_numeric(incidence["resistant_infection_n"], errors="coerce")
        pop = pd.to_numeric(incidence["population_n"], errors="coerce")
        too_many = _idx(incidence, cases.notna() & pop.notna() & (cases > pop))
        if too_many:
            f.add(
                "incidence",
                "ERROR",
                "cases_exceed_population",
                "resistant_infection_n cannot exceed population_n.",
                too_many,
            )
    if "population_n" in incidence.columns:
        population = pd.to_numeric(incidence["population_n"], errors="coerce")
        zero_population = _idx(incidence, population.notna() & (population == 0))
        if zero_population:
            f.add(
                "incidence",
                "ERROR",
                "zero_population_denominator",
                "Population denominator must be greater than zero.",
                zero_population,
            )
    for col in ("source_name", "case_definition"):
        if col in incidence.columns:
            missing = _idx(incidence, incidence[col].isna())
            if missing:
                f.add(
                    "incidence",
                    "WARNING",
                    f"missing_{col}",
                    f"Incidence rows without {col} have weak provenance or an ambiguous numerator.",
                    missing,
                )
    _duplicate_key_check(
        f,
        "incidence",
        incidence,
        [
            "year",
            "country",
            "source_sector",
            "organism",
            "antimicrobial",
            "source_name",
            "case_definition",
        ],
        "duplicate_incidence_key",
        "incidence",
    )


def _audit_consumption(f: _Findings, consumption: pd.DataFrame) -> None:
    if consumption.empty:
        return
    numeric_values: dict[str, pd.Series] = {}
    for col in ("ddd_total", "population_n", "period_days"):
        if col not in consumption.columns:
            continue
        values = pd.to_numeric(consumption[col], errors="coerce")
        numeric_values[col] = values
        missing = _idx(consumption, values.isna())
        if missing:
            f.add(
                "consumption",
                "ERROR",
                f"missing_{col}",
                f"{col} is required to calculate DDD/1,000 inhabitants/day.",
                missing,
            )
        negative = _idx(consumption, values.notna() & (values < 0))
        if negative:
            f.add("consumption", "ERROR", f"negative_{col}", f"{col} cannot be negative.", negative)
        infinite = _idx(consumption, values.notna() & ~np.isfinite(values))
        if infinite:
            f.add("consumption", "ERROR", f"non_finite_{col}", f"{col} must be finite.", infinite)
    for col in ("population_n", "period_days"):
        values = numeric_values.get(col)
        if values is not None:
            fractional = _idx(consumption, values.notna() & (values % 1 != 0))
            if fractional:
                f.add(
                    "consumption",
                    "ERROR",
                    f"non_integer_{col}",
                    f"{col} must be a whole number.",
                    fractional,
                )
            zero = _idx(consumption, values.notna() & (values == 0))
            if zero:
                f.add(
                    "consumption", "ERROR", f"zero_{col}", f"{col} must be greater than zero.", zero
                )
    _audit_years(f, "consumption", consumption)
    days = numeric_values.get("period_days")
    if days is not None:
        long_period = _idx(consumption, days.notna() & (days > 366))
        if long_period:
            f.add(
                "consumption",
                "ERROR",
                "implausible_period_days",
                "period_days exceeds 366; one row must describe at most one year.",
                long_period,
            )
    if "care_sector" in consumption.columns:
        sector = consumption["care_sector"].dropna().astype(str).str.lower()
        unknown_sector = sector[~sector.isin(VALID_CARE_SECTORS)].unique().tolist()
        if unknown_sector:
            f.add(
                "consumption",
                "WARNING",
                "nonstandard_care_sector",
                "Unrecognised care-sector labels are retained; document the scope of the consumption source.",
                unknown_sector,
            )
    if "data_source_type" in consumption.columns:
        kinds = consumption["data_source_type"].dropna().astype(str).str.lower()
        unknown = kinds[~kinds.isin(VALID_CONSUMPTION_SOURCES)].unique().tolist()
        if unknown:
            f.add(
                "consumption",
                "WARNING",
                "nonstandard_consumption_source_type",
                "Unrecognised source types are retained; document how consumption was measured.",
                unknown,
            )
    if "source_name" in consumption.columns:
        missing = _idx(consumption, consumption["source_name"].isna())
        if missing:
            f.add(
                "consumption",
                "WARNING",
                "missing_source_name",
                "Consumption rows without a source name have weak provenance.",
                missing,
            )
    _duplicate_key_check(
        f,
        "consumption",
        consumption,
        [
            "year",
            "country",
            "care_sector",
            "antimicrobial",
            "atc_code",
            "source_name",
            "data_source_type",
            "period_days",
        ],
        "duplicate_consumption_key",
        "consumption",
    )
