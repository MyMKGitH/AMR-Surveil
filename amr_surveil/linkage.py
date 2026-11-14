"""Exact-isolate linkage utilities. No fuzzy matching or ecological joins are permitted."""

from __future__ import annotations

import json
from collections import defaultdict

import pandas as pd

from amr_surveil.manifest import unique_manifest

_unique_manifest = unique_manifest


def build_matched_isolates(
    samples: pd.DataFrame,
    amr_calls: pd.DataFrame,
    ast: pd.DataFrame,
) -> pd.DataFrame:
    """Summarise isolates having both determinant calls and AST on exact sample_id.

    This creates a joined evidence inventory, not genotype-to-phenotype prediction.
    It does not assume every detected gene is relevant to every antimicrobial.
    """
    columns = [
        "sample_id",
        "organism",
        "country",
        "source_sector",
        "collection_date",
        "n_detected_determinants",
        "detected_determinants",
        "n_ast_records",
        "ast_evidence",
        "breakpoint_context_complete",
        "interpretation_boundary",
    ]
    if samples.empty or amr_calls.empty or ast.empty:
        return pd.DataFrame(columns=columns)

    manifest = _unique_manifest(samples)
    calls = amr_calls.dropna(subset=["sample_id", "gene"]).copy()
    ast_valid = ast.copy()
    ast_valid["result"] = ast_valid["result"].fillna("").str.upper()
    ast_valid = ast_valid[ast_valid["result"].isin(["S", "I", "R"])].copy()
    if manifest.empty or calls.empty or ast_valid.empty:
        return pd.DataFrame(columns=columns)

    known = set(manifest["sample_id"].astype(str))
    calls = calls[calls["sample_id"].astype(str).isin(known)].copy()
    ast_valid = ast_valid[ast_valid["sample_id"].astype(str).isin(known)].copy()
    if calls.empty or ast_valid.empty:
        return pd.DataFrame(columns=columns)

    # An explicit organism disagreement in either isolate-level stream excludes the sample.
    mismatch_ids: set[str] = set()
    if "organism" in manifest.columns:
        sample_org = manifest[["sample_id", "organism"]].dropna().drop_duplicates()
        for evidence in (calls, ast_valid):
            if "organism" not in evidence.columns:
                continue
            evidence_org = evidence[["sample_id", "organism"]].dropna().drop_duplicates()
            compare = evidence_org.merge(sample_org, on="sample_id", suffixes=("_input", "_sample"))
            input_org = compare["organism_input"].astype("string").str.strip().str.casefold()
            manifest_org = compare["organism_sample"].astype("string").str.strip().str.casefold()
            mismatches = compare.loc[input_org != manifest_org, "sample_id"].astype(str)
            mismatch_ids.update(mismatches.tolist())
    calls = calls[~calls["sample_id"].astype(str).isin(mismatch_ids)].copy()
    ast_valid = ast_valid[~ast_valid["sample_id"].astype(str).isin(mismatch_ids)].copy()
    if calls.empty or ast_valid.empty:
        return pd.DataFrame(columns=columns)

    determinant_sets: dict[str, set[str]] = defaultdict(set)
    for sample_id, gene in calls[["sample_id", "gene"]].itertuples(index=False, name=None):
        if pd.notna(sample_id) and pd.notna(gene) and str(gene).strip():
            determinant_sets[str(sample_id)].add(str(gene))
    detections = pd.DataFrame(
        [
            {
                "sample_id": sample_id,
                "n_detected_determinants": len(genes),
                "detected_determinants": json.dumps(sorted(genes), ensure_ascii=False),
            }
            for sample_id, genes in sorted(determinant_sets.items())
        ]
    )

    ast_fields = [
        col
        for col in [
            "organism",
            "antimicrobial",
            "result",
            "mic",
            "mic_unit",
            "method",
            "breakpoint_standard",
            "breakpoint_version",
        ]
        if col in ast_valid.columns
    ]
    context_fields = [
        c for c in ["method", "breakpoint_standard", "breakpoint_version"] if c in ast_valid.columns
    ]
    ast_payloads: dict[str, list[dict[str, object]]] = defaultdict(list)
    ast_context: dict[str, bool] = defaultdict(lambda: bool(context_fields))
    ast_columns = ["sample_id", *ast_fields]
    for values in ast_valid[ast_columns].itertuples(index=False, name=None):
        sample_id = str(values[0])
        record: dict[str, object] = {}
        for col, value in zip(ast_fields, values[1:], strict=True):
            record[col] = "" if pd.isna(value) else str(value)
        ast_payloads[sample_id].append(record)
        if any(not record.get(col) for col in context_fields):
            ast_context[sample_id] = False
    ast_records = [
        {
            "sample_id": sample_id,
            "n_ast_records": len(payload),
            "ast_evidence": json.dumps(payload, ensure_ascii=False, sort_keys=True),
            "breakpoint_context_complete": ast_context[sample_id],
        }
        for sample_id, payload in sorted(ast_payloads.items())
    ]
    phenotypes = pd.DataFrame(ast_records)
    paired = detections.merge(phenotypes, on="sample_id", how="inner", validate="one_to_one")
    meta_cols = [
        col
        for col in ["sample_id", "organism", "country", "source_sector", "collection_date"]
        if col in manifest.columns
    ]
    paired = paired.merge(manifest[meta_cols], on="sample_id", how="inner", validate="one_to_one")
    for col in columns:
        if col not in paired.columns:
            paired[col] = pd.NA
    paired["interpretation_boundary"] = (
        "Exact sample-ID co-occurrence only; gene-to-drug relevance and clinical susceptibility are not inferred."
    )
    return paired[columns].sort_values("sample_id").reset_index(drop=True)


def build_prediction_context(samples: pd.DataFrame, predictions: pd.DataFrame) -> pd.DataFrame:
    """Attach metadata to model outputs by exact sample_id, keeping predictions separate from AST."""
    if predictions.empty or samples.empty:
        return pd.DataFrame(
            columns=list(predictions.columns) + ["organism", "source_sector", "collection_date"]
        )
    manifest = _unique_manifest(samples)
    metadata_cols = [
        col
        for col in ["sample_id", "organism", "source_sector", "collection_date", "country"]
        if col in manifest.columns
    ]
    pred = predictions[
        predictions["sample_id"].astype(str).isin(set(manifest["sample_id"].astype(str)))
    ].copy()
    merged = pred.merge(
        manifest[metadata_cols],
        on="sample_id",
        how="inner",
        suffixes=("", "_sample"),
        validate="many_to_one",
    )
    if "organism_sample" in merged.columns:
        if "organism" in merged.columns:
            model_org = merged["organism"].astype("string").str.strip().str.casefold()
            sample_org = merged["organism_sample"].astype("string").str.strip().str.casefold()
            mismatch = model_org.notna() & sample_org.notna() & (model_org != sample_org)
            merged = merged[~mismatch].copy()
            merged["organism"] = merged["organism"].fillna(merged["organism_sample"])
            merged = merged.drop(columns=["organism_sample"])
        else:
            merged = merged.rename(columns={"organism_sample": "organism"})
    return merged


def gene_drug_stratification(
    calls: pd.DataFrame,
    ast: pd.DataFrame,
    mapping: pd.DataFrame,
    samples: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Descriptively stratify AST among gene-positive isolates using a curated map.

    Gene, organism, exact sample ID, antimicrobial, method and breakpoint context are
    required to form a descriptive stratum. This is not diagnostic accuracy analysis.
    """
    columns = [
        "organism",
        "gene",
        "antimicrobial",
        "mapping_source",
        "method",
        "breakpoint_standard",
        "breakpoint_version",
        "n_gene_positive_with_ast",
        "n_R",
        "n_I",
        "n_S",
        "R_pct_among_gene_positive_with_ast",
        "interpretation_boundary",
    ]
    if calls.empty or ast.empty or mapping.empty or samples is None or samples.empty:
        return pd.DataFrame(columns=columns)
    required = {"organism", "gene", "antimicrobial", "mapping_source"}
    if not required.issubset(mapping.columns):
        return pd.DataFrame(columns=columns)

    manifest = _unique_manifest(samples)
    if manifest.empty or "organism" not in manifest.columns:
        return pd.DataFrame(columns=columns)
    manifest = manifest[["sample_id", "organism"]].copy()
    manifest["organism_manifest"] = manifest["organism"].astype("string").str.strip()
    manifest = manifest.drop(columns="organism")

    map_table = (
        mapping.dropna(subset=["organism", "gene", "antimicrobial", "mapping_source"])
        .copy()
        .drop_duplicates()
    )
    for col in ["organism", "gene", "antimicrobial", "mapping_source"]:
        map_table[col] = map_table[col].astype("string").str.strip()

    # Enrich missing organism labels from the unique manifest, but never overwrite
    # a conflicting organism label supplied by an upstream detector.
    call_table = calls.dropna(subset=["sample_id", "gene"]).copy()
    call_keep = [c for c in ["sample_id", "gene", "organism"] if c in call_table.columns]
    call_table = call_table[call_keep].drop_duplicates()
    call_table = call_table.merge(manifest, on="sample_id", how="inner", validate="many_to_one")
    if "organism" not in call_table.columns:
        call_table["organism"] = call_table["organism_manifest"]
    else:
        supplied = call_table["organism"].astype("string").str.strip()
        expected = call_table["organism_manifest"].astype("string").str.strip()
        mismatch = (
            supplied.notna()
            & expected.notna()
            & (supplied.str.casefold() != expected.str.casefold())
        )
        call_table = call_table[~mismatch].copy()
        call_table["organism"] = call_table["organism"].fillna(call_table["organism_manifest"])
    call_table = call_table.drop(columns="organism_manifest")
    call_table = call_table.dropna(subset=["organism", "gene"]).drop_duplicates()

    # AST organism is also enriched from the manifest if absent; a conflict is excluded.
    ast_table = ast.dropna(subset=["sample_id", "antimicrobial", "result"]).copy()
    ast_table["result"] = ast_table["result"].astype("string").str.strip().str.upper()
    ast_table = ast_table[ast_table["result"].isin(["S", "I", "R"])].copy()
    ast_table = ast_table.merge(manifest, on="sample_id", how="inner", validate="many_to_one")
    if "organism" not in ast_table.columns:
        ast_table["organism"] = ast_table["organism_manifest"]
    else:
        supplied = ast_table["organism"].astype("string").str.strip()
        expected = ast_table["organism_manifest"].astype("string").str.strip()
        mismatch = (
            supplied.notna()
            & expected.notna()
            & (supplied.str.casefold() != expected.str.casefold())
        )
        ast_table = ast_table[~mismatch].copy()
        ast_table["organism"] = ast_table["organism"].fillna(ast_table["organism_manifest"])
    ast_table = ast_table.drop(columns="organism_manifest")
    ast_table = ast_table.dropna(subset=["organism", "antimicrobial"])

    joined = call_table.merge(map_table, on=["gene", "organism"], how="inner").merge(
        ast_table,
        on=["sample_id", "organism", "antimicrobial"],
        how="inner",
        suffixes=("_map", "_ast"),
    )
    if joined.empty:
        return pd.DataFrame(columns=columns)

    group_cols = [
        "organism",
        "gene",
        "antimicrobial",
        "mapping_source",
        "method",
        "breakpoint_standard",
        "breakpoint_version",
    ]
    for col in group_cols:
        if col not in joined.columns:
            joined[col] = "Unspecified"
        joined[col] = joined[col].fillna("Unspecified")

    rows: list[dict[str, object]] = []
    for key, group in joined.groupby(group_cols, dropna=False, sort=True):
        # Conflicting observations in the same isolate/stratum are excluded instead
        # of selecting an arbitrary category in exploratory mode.
        conflicting = group.groupby("sample_id")["result"].nunique()
        ambiguous_ids = set(conflicting[conflicting > 1].index.astype(str))
        group = group[~group["sample_id"].astype(str).isin(ambiguous_ids)]
        group = group.drop_duplicates("sample_id")
        n = int(group["sample_id"].nunique())
        if not n:
            continue
        counts = group["result"].value_counts()
        row = dict(zip(group_cols, key, strict=True))
        row.update(
            {
                "n_gene_positive_with_ast": n,
                "n_R": int(counts.get("R", 0)),
                "n_I": int(counts.get("I", 0)),
                "n_S": int(counts.get("S", 0)),
                "R_pct_among_gene_positive_with_ast": round(100 * int(counts.get("R", 0)) / n, 4),
                "interpretation_boundary": (
                    "Descriptive among gene-positive isolates only; not sensitivity, specificity, "
                    "PPV or a clinical prediction."
                ),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)
