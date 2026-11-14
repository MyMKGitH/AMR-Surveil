"""Genotype-phenotype readiness gate and descriptive concordance.

Go / no-go rule: proceed with a
genotype-phenotype analysis only if WGS and AST can be verified to belong to the same isolate,
the AST method and interpretive standard are identified, and the sample-selection process is
understood. :func:`assess_readiness` encodes that rule; :func:`genotype_phenotype_concordance`
refuses to run unless it returns ``GO``.

Scientific boundaries enforced here
-----------------------------------
* A "negative" genotype is defined **only** for isolates declared as screened (screening table)
  with a single, recorded method/database/version. A missing call row for an unscreened isolate is
  never a negative.
* Gene-to-drug relationships come exclusively from a user-supplied, cited map; none is built in.
* Output is descriptive concordance between a mapped-determinant screen and the recorded
  phenotype. It is not diagnostic accuracy, a clinical prediction or a causal attribution.
* ``I`` (EUCAST: susceptible, increased exposure) is counted as *not resistant* and reported
  separately (``n_I_in_notR``).
"""

from __future__ import annotations

import json
from typing import Any

import pandas as pd

from amr_surveil.epidemiology import wilson_interval
from amr_surveil.exclusions import ExclusionLog
from amr_surveil.manifest import unique_manifest
from amr_surveil.normalize import label_key

CONCORDANCE_COLUMNS = [
    "organism",
    "antimicrobial",
    "method",
    "breakpoint_standard",
    "breakpoint_version",
    "screen_database",
    "screen_database_version",
    "mapping_sources",
    "n_isolates",
    "n_gene_pos_R",
    "n_gene_pos_notR",
    "n_gene_neg_R",
    "n_gene_neg_notR",
    "n_I_in_notR",
    "n_excluded_uncertain_genotype",
    "sensitivity_pct",
    "sensitivity_low_pct",
    "sensitivity_high_pct",
    "specificity_pct",
    "specificity_low_pct",
    "specificity_high_pct",
    "ppv_pct",
    "ppv_low_pct",
    "ppv_high_pct",
    "npv_pct",
    "npv_low_pct",
    "npv_high_pct",
    "metrics_suppressed",
    "interpretation_boundary",
]
DISCORDANT_COLUMNS = [
    "sample_id",
    "organism",
    "antimicrobial",
    "discordance_type",
    "mapped_genes_detected",
    "result",
    "mic",
    "method",
    "breakpoint_standard",
    "breakpoint_version",
    "screen_database",
    "screen_database_version",
]
BOUNDARY = (
    "Descriptive concordance of a mapped-determinant screen with the recorded phenotype among screened "
    "isolates; not diagnostic accuracy, causal attribution or a clinical prediction. A negative screen is "
    "not evidence of susceptibility. I is counted as not-R."
)


def _valid_ast(ast: pd.DataFrame) -> pd.DataFrame:
    if ast.empty:
        return ast.iloc[0:0].copy()
    out = ast.copy()
    out["result"] = out["result"].astype("string").str.strip().str.upper()
    return out[out["result"].isin(["S", "I", "R"])].copy()


def _single_screening(
    screening: pd.DataFrame, manifest_ids: set[str], exclusions: ExclusionLog | None
):
    """Screened samples with exactly one distinct method/database/version."""
    cols = ["method", "database", "database_version"]
    if screening.empty or "sample_id" not in screening.columns:
        return pd.DataFrame(
            columns=["sample_id", "screen_method", "screen_database", "screen_database_version"]
        )
    scr = screening.dropna(subset=["sample_id"]).copy()
    scr["sample_id"] = scr["sample_id"].astype(str)
    for c in cols:
        if c not in scr.columns:
            scr[c] = pd.NA
    scr = scr[scr["sample_id"].isin(manifest_ids)]
    distinct = scr[["sample_id", *cols]].drop_duplicates()
    counts = distinct["sample_id"].value_counts()
    multi = set(counts[counts > 1].index)
    if exclusions is not None:
        for sid in sorted(multi):
            exclusions.add(
                "concordance",
                "screening",
                sid,
                "multiple_screening_records",
                "more than one distinct screening method/database/version",
            )
    single = distinct[~distinct["sample_id"].isin(multi)].rename(
        columns={
            "method": "screen_method",
            "database": "screen_database",
            "database_version": "screen_database_version",
        }
    )
    return single.reset_index(drop=True)


def assess_readiness(
    samples: pd.DataFrame,
    amr_calls: pd.DataFrame,
    ast: pd.DataFrame,
    screening: pd.DataFrame,
    gene_map: pd.DataFrame,
    study: dict[str, Any],
) -> dict[str, Any]:
    """Evaluate the go/no-go criteria. Returns decision, criteria and blocking reasons."""
    manifest = unique_manifest(samples)
    manifest_ids = set(manifest["sample_id"].astype(str)) if not manifest.empty else set()
    single = _single_screening(screening, manifest_ids, None)
    valid = _valid_ast(ast)
    paired_ast = (
        valid[valid["sample_id"].astype(str).isin(set(single["sample_id"]))]
        if not valid.empty
        else valid
    )
    paired_ids = set(paired_ast["sample_id"].astype(str)) if not paired_ast.empty else set()

    def all_present(frame: pd.DataFrame, cols: list[str]) -> bool:
        return (not frame.empty) and all(
            c in frame.columns and frame[c].notna().all() for c in cols
        )

    scr_paired = single[single["sample_id"].isin(paired_ids)]
    org_conflict = 0
    if not paired_ast.empty and "organism" in paired_ast.columns and "organism" in manifest.columns:
        comp = (
            paired_ast[["sample_id", "organism"]]
            .dropna()
            .drop_duplicates()
            .merge(
                manifest[["sample_id", "organism"]].dropna(), on="sample_id", suffixes=("_a", "_m")
            )
        )
        org_conflict = int(
            (comp["organism_a"].map(label_key) != comp["organism_m"].map(label_key)).sum()
        )

    map_ok = (
        (not gene_map.empty) and gene_map["mapping_source"].notna().all()
        if "mapping_source" in gene_map.columns
        else False
    )
    attested = bool(study.get("sample_selection_documented")) and bool(
        str(study.get("selection_description", "")).strip()
    )

    criteria = [
        (
            "paired_isolate_level",
            "WGS screening and AST share the same isolate (exact sample_id)",
            len(paired_ids) > 0,
            f"{len(paired_ids)} isolate(s) have a single screening record and at least one valid S/I/R result",
        ),
        (
            "ast_method_recorded",
            "AST method recorded for every paired AST record",
            all_present(paired_ast, ["method"]),
            "method must be populated",
        ),
        (
            "interpretive_standard_recorded",
            "Breakpoint standard and version recorded for every paired AST record",
            all_present(paired_ast, ["breakpoint_standard", "breakpoint_version"]),
            "breakpoint_standard and breakpoint_version must be populated",
        ),
        (
            "screening_provenance_recorded",
            "Screening tool, database and database version recorded for every paired isolate",
            all_present(
                scr_paired, ["screen_method", "screen_database", "screen_database_version"]
            ),
            "screening table needs method, database and database_version",
        ),
        (
            "organism_consistent",
            "No organism conflicts between manifest and AST among paired isolates",
            org_conflict == 0,
            f"{org_conflict} conflicting isolate/organism label(s)",
        ),
        (
            "cited_gene_drug_map",
            "A gene-to-drug map with a source citation for every row was supplied",
            bool(map_ok),
            "supply --gene-drug-map with mapping_source for every row",
        ),
        (
            "sample_selection_documented",
            "Sample-selection process documented by the analyst (cannot be verified by software)",
            attested,
            "set study.sample_selection_documented: true and describe it in study.selection_description",
        ),
    ]
    rows = [{"id": i, "description": d, "met": bool(m), "detail": t} for i, d, m, t in criteria]
    blocking = [r["description"] + " — " + r["detail"] for r in rows if not r["met"]]
    return {
        "decision": "GO" if not blocking else "NO-GO",
        "criteria": rows,
        "blocking_reasons": blocking,
        "n_paired_isolates": len(paired_ids),
    }


def _metric(num: int, den: int, confidence_level: float) -> tuple[float, float, float]:
    if den <= 0:
        return (float("nan"),) * 3
    low, high = wilson_interval(num, den, confidence_level)
    return round(100 * num / den, 4), round(low, 4), round(high, 4)


def _match_positives(calls: pd.DataFrame, gene_map: pd.DataFrame) -> pd.DataFrame:
    """(sample_id, drug_key, org_key, gene) for every call matching a map row."""
    pieces = []
    cm = gene_map.copy()
    if "match_type" not in cm.columns:
        cm["match_type"] = "exact"
    cm["match_type"] = cm["match_type"].fillna("exact").astype(str).str.strip().str.lower()
    cm["org_key"] = cm["organism"].map(label_key)
    cm["drug_key"] = cm["antimicrobial"].map(label_key)
    cm["map_gene"] = cm["gene"].astype(str).str.strip()
    exact = cm[cm["match_type"] != "prefix"]
    if not exact.empty:
        pieces.append(
            calls.merge(
                exact[["org_key", "map_gene", "drug_key"]],
                left_on=["org_key", "gene"],
                right_on=["org_key", "map_gene"],
                how="inner",
            )[["sample_id", "org_key", "drug_key", "gene"]]
        )
    for row in cm[cm["match_type"] == "prefix"].itertuples(index=False):
        sub = calls[calls["org_key"] == row.org_key]
        hit = sub[(sub["gene"] == row.map_gene) | sub["gene"].str.startswith(row.map_gene + "-")]
        if not hit.empty:
            piece = hit[["sample_id", "org_key", "gene"]].copy()
            piece["drug_key"] = row.drug_key
            pieces.append(piece[["sample_id", "org_key", "drug_key", "gene"]])
    if not pieces:
        return pd.DataFrame(columns=["sample_id", "org_key", "drug_key", "gene"])
    return pd.concat(pieces, ignore_index=True).drop_duplicates()


def genotype_phenotype_concordance(
    samples: pd.DataFrame,
    amr_calls: pd.DataFrame,
    ast: pd.DataFrame,
    screening: pd.DataFrame,
    gene_map: pd.DataFrame,
    *,
    uncertain_calls: pd.DataFrame | None = None,
    min_isolates: int = 10,
    confidence_level: float = 0.95,
    exclusions: ExclusionLog | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """2x2 concordance per organism/drug/method/breakpoint/screening stratum, plus discordants.

    ``amr_calls`` must already be restricted to the evidence tiers the analysis accepts. Callers
    should check :func:`assess_readiness` first; the pipeline does so and skips this function on NO-GO.

    ``uncertain_calls`` are calls that matched a mapped gene but were set aside (weak evidence tier,
    partial or internal-stop hit). An isolate whose *only* mapped hit for a drug is uncertain is
    neither genotype-positive nor genotype-negative for that drug, so it is left out of the 2x2 table,
    logged, and counted in ``n_excluded_uncertain_genotype`` instead of being counted as a negative.
    """
    empty = (pd.DataFrame(columns=CONCORDANCE_COLUMNS), pd.DataFrame(columns=DISCORDANT_COLUMNS))
    manifest = unique_manifest(samples)
    if manifest.empty or "organism" not in manifest.columns or gene_map.empty or ast.empty:
        return empty
    manifest = manifest[["sample_id", "organism"]].copy()
    manifest["sample_id"] = manifest["sample_id"].astype(str)
    manifest["org_key"] = manifest["organism"].map(label_key)
    manifest = manifest.rename(columns={"organism": "organism_manifest"})
    single = _single_screening(screening, set(manifest["sample_id"]), exclusions)
    if single.empty:
        return empty
    screened = set(single["sample_id"])

    # --- phenotype side
    valid = _valid_ast(ast)
    valid = valid[valid["sample_id"].astype(str).isin(screened)].copy()
    valid["sample_id"] = valid["sample_id"].astype(str)
    valid = valid.merge(manifest, on="sample_id", how="inner", validate="many_to_one")
    if "organism" in valid.columns:
        ast_org = valid["organism"].map(label_key)
        conflict = ast_org.notna() & (ast_org != valid["org_key"])
        if exclusions is not None:
            for sid in valid.loc[conflict, "sample_id"].unique():
                exclusions.add(
                    "concordance",
                    "ast",
                    sid,
                    "organism_mismatch",
                    "AST organism conflicts with manifest",
                )
        valid = valid[~conflict].drop(columns="organism")
    valid = valid.rename(columns={"organism_manifest": "organism"})
    for col in ("method", "breakpoint_standard", "breakpoint_version", "mic"):
        if col not in valid.columns:
            valid[col] = pd.NA
    valid["drug_key"] = valid["antimicrobial"].map(label_key)
    key = ["sample_id", "drug_key", "method", "breakpoint_standard", "breakpoint_version"]
    nun = valid.groupby(key, dropna=False)["result"].transform("nunique")
    if exclusions is not None:
        for sid, drug in (
            valid.loc[nun > 1, ["sample_id", "antimicrobial"]]
            .drop_duplicates()
            .itertuples(index=False)
        ):
            exclusions.add(
                "concordance",
                "ast",
                f"{sid}|{drug}",
                "conflicting_duplicate_ast",
                "conflicting S/I/R",
            )
    valid = valid[nun == 1].drop_duplicates(key + ["result"])

    # --- mapped drug universe: only organism/drug pairs present in the cited map
    gm = gene_map.dropna(subset=["organism", "gene", "antimicrobial", "mapping_source"]).copy()
    gm["org_key"] = gm["organism"].map(label_key)
    gm["drug_key"] = gm["antimicrobial"].map(label_key)
    sources = (
        gm.groupby(["org_key", "drug_key"])["mapping_source"]
        .apply(lambda s: "; ".join(sorted(set(map(str, s)))))
        .rename("mapping_sources")
        .reset_index()
    )
    valid = valid.merge(sources, on=["org_key", "drug_key"], how="inner")
    if valid.empty:
        return empty

    # --- genotype side
    calls = amr_calls.dropna(subset=["sample_id", "gene"]).copy()
    calls["sample_id"] = calls["sample_id"].astype(str)
    calls = calls[calls["sample_id"].isin(screened)].merge(manifest, on="sample_id", how="inner")
    if "organism" in calls.columns:
        c_org = calls["organism"].map(label_key)
        calls = calls[~(c_org.notna() & (c_org != calls["org_key"]))].drop(columns="organism")
    calls["gene"] = calls["gene"].astype(str).str.strip()
    positives = _match_positives(calls, gm)
    pos_genes = (
        positives.groupby(["sample_id", "drug_key"])["gene"]
        .apply(lambda s: json.dumps(sorted(set(s))))
        .rename("mapped_genes_detected")
        .reset_index()
    )
    valid = valid.merge(pos_genes, on=["sample_id", "drug_key"], how="left")
    valid["gene_pos"] = valid["mapped_genes_detected"].notna()
    uncertain_pairs: set[tuple[str, str]] = set()
    if uncertain_calls is not None and not uncertain_calls.empty:
        unc = uncertain_calls.dropna(subset=["sample_id", "gene"]).copy()
        unc["sample_id"] = unc["sample_id"].astype(str)
        unc = unc[unc["sample_id"].isin(screened)].merge(manifest, on="sample_id", how="inner")
        if "organism" in unc.columns:
            u_org = unc["organism"].map(label_key)
            unc = unc[~(u_org.notna() & (u_org != unc["org_key"]))].drop(columns="organism")
        unc["gene"] = unc["gene"].astype(str).str.strip()
        hits = _match_positives(unc, gm)
        uncertain_pairs = set(zip(hits["sample_id"], hits["drug_key"], strict=True))
    valid["uncertain"] = [
        (not pos) and (sid, drug) in uncertain_pairs
        for sid, drug, pos in zip(
            valid["sample_id"], valid["drug_key"], valid["gene_pos"], strict=True
        )
    ]
    if exclusions is not None:
        for sid, drug in valid.loc[valid["uncertain"], ["sample_id", "antimicrobial"]].itertuples(
            index=False
        ):
            exclusions.add(
                "concordance",
                "ast",
                f"{sid}|{drug}",
                "uncertain_genotype_low_evidence_hit",
                "only a weak-tier or non-functional hit for a mapped gene; not counted as negative",
            )
    valid["is_R"] = valid["result"] == "R"
    valid = valid.merge(
        single[["sample_id", "screen_database", "screen_database_version"]],
        on="sample_id",
        how="left",
    )

    group_cols = [
        "organism",
        "antimicrobial",
        "method",
        "breakpoint_standard",
        "breakpoint_version",
        "screen_database",
        "screen_database_version",
        "mapping_sources",
    ]
    for col in group_cols:
        valid[col] = valid[col].fillna("Unspecified")
    rows = []
    for gkey, g_all in valid.groupby(group_cols, sort=True, dropna=False):
        g = g_all[~g_all["uncertain"]]
        tp = int((g.gene_pos & g.is_R).sum())
        fp = int((g.gene_pos & ~g.is_R).sum())
        fn = int((~g.gene_pos & g.is_R).sum())
        tn = int((~g.gene_pos & ~g.is_R).sum())
        n = tp + fp + fn + tn
        suppressed = n < min_isolates
        row = dict(zip(group_cols, gkey, strict=True))
        row.update(
            n_isolates=n,
            n_gene_pos_R=tp,
            n_gene_pos_notR=fp,
            n_gene_neg_R=fn,
            n_gene_neg_notR=tn,
            n_I_in_notR=int((~g.is_R & (g.result == "I")).sum()),
            n_excluded_uncertain_genotype=int(g_all["uncertain"].sum()),
            metrics_suppressed=suppressed,
            interpretation_boundary=BOUNDARY,
        )
        for name, (num, den) in {
            "sensitivity": (tp, tp + fn),
            "specificity": (tn, tn + fp),
            "ppv": (tp, tp + fp),
            "npv": (tn, tn + fn),
        }.items():
            pct, low, high = (
                (float("nan"),) * 3 if suppressed else _metric(num, den, confidence_level)
            )
            row[f"{name}_pct"], row[f"{name}_low_pct"], row[f"{name}_high_pct"] = pct, low, high
        rows.append(row)
    summary = pd.DataFrame(rows, columns=CONCORDANCE_COLUMNS)

    valid = valid[~valid["uncertain"]]
    disc = valid[(valid.gene_pos & ~valid.is_R) | (~valid.gene_pos & valid.is_R)].copy()
    disc["discordance_type"] = (
        disc.apply(
            lambda r: (
                "mapped_gene_detected_but_phenotype_not_R"
                if r["gene_pos"]
                else "phenotype_R_without_mapped_gene"
            ),
            axis=1,
        )
        if not disc.empty
        else pd.Series(dtype=str)
    )
    disc["mapped_genes_detected"] = disc["mapped_genes_detected"].fillna("[]")
    for col in DISCORDANT_COLUMNS:
        if col not in disc.columns:
            disc[col] = pd.NA
    return summary, disc[DISCORDANT_COLUMNS].sort_values(
        ["sample_id", "antimicrobial"]
    ).reset_index(drop=True)
