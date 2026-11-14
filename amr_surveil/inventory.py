"""Data inventory and AST missingness summary."""

from __future__ import annotations

import pandas as pd

from amr_surveil.manifest import ambiguous_sample_ids
from amr_surveil.normalize import parse_iso_dates

INVENTORY_COLUMNS = [
    "sample_id",
    "sequencing_accession",
    "source_study",
    "data_license",
    "organism",
    "collection_year",
    "country",
    "region",
    "source_sector",
    "specimen",
    "sequencing_type",
    "amr_screened",
    "n_amr_calls",
    "n_plasmid_calls",
    "n_ast_records",
    "n_predictions",
    "has_genotype_and_ast",
    "ambiguous_sample_id",
]
MISSINGNESS_COLUMNS = [
    "country",
    "organism",
    "antimicrobial",
    "year",
    "testing_lab",
    "n_records",
    "n_valid_result",
    "n_missing_or_invalid_result",
    "pct_missing_or_invalid",
]


def _count(frame: pd.DataFrame) -> pd.Series:
    if frame.empty or "sample_id" not in frame.columns:
        return pd.Series(dtype="int64")
    return frame["sample_id"].dropna().astype(str).value_counts()


def build_data_inventory(
    samples: pd.DataFrame,
    amr_calls: pd.DataFrame,
    screening: pd.DataFrame,
    ast: pd.DataFrame,
    plasmids: pd.DataFrame,
    predictions: pd.DataFrame,
) -> pd.DataFrame:
    """One row per manifest row: provenance fields plus which evidence streams exist."""
    if samples.empty:
        return pd.DataFrame(columns=INVENTORY_COLUMNS)
    inv = samples.copy()
    sid = inv["sample_id"].astype("string")
    calls, asts, plas, preds = _count(amr_calls), _count(ast), _count(plasmids), _count(predictions)
    screened = set(screening["sample_id"].dropna().astype(str)) if not screening.empty else set()
    valid_ast = set()
    if not ast.empty:
        res = ast["result"].astype("string").str.strip().str.upper()
        valid_ast = set(ast.loc[res.isin(["S", "I", "R"]), "sample_id"].dropna().astype(str))
    inv["collection_year"] = (
        parse_iso_dates(inv["collection_date"])["year"]
        if "collection_date" in inv.columns
        else pd.NA
    )
    inv["amr_screened"] = sid.isin(screened)
    inv["n_amr_calls"] = sid.map(calls).fillna(0).astype(int)
    inv["n_plasmid_calls"] = sid.map(plas).fillna(0).astype(int)
    inv["n_ast_records"] = sid.map(asts).fillna(0).astype(int)
    inv["n_predictions"] = sid.map(preds).fillna(0).astype(int)
    inv["has_genotype_and_ast"] = inv["amr_screened"] & sid.isin(valid_ast)
    inv["ambiguous_sample_id"] = sid.isin(ambiguous_sample_ids(samples))
    for col in INVENTORY_COLUMNS:
        if col not in inv.columns:
            inv[col] = pd.NA
    return inv[INVENTORY_COLUMNS]


def missingness_summary(samples: pd.DataFrame, ast: pd.DataFrame) -> pd.DataFrame:
    """AST result completeness by country, organism, antimicrobial, year and laboratory."""
    if ast.empty:
        return pd.DataFrame(columns=MISSINGNESS_COLUMNS)
    work = ast.copy()
    meta = samples.drop_duplicates("sample_id", keep=False)
    cols = [c for c in ["sample_id", "country", "organism", "collection_date"] if c in meta.columns]
    work = work.merge(meta[cols], on="sample_id", how="left", suffixes=("", "_manifest"))
    org = work["organism"] if "organism" in work.columns else pd.Series(pd.NA, index=work.index)
    if "organism_manifest" in work.columns:
        org = org.fillna(work["organism_manifest"])
    work["organism"] = org
    work["country"] = work["country"].fillna("Unknown") if "country" in work.columns else "Unknown"
    work["year"] = (
        parse_iso_dates(work["collection_date"])["year"].astype("string")
        if "collection_date" in work.columns
        else pd.NA
    )
    work["antimicrobial"] = work["antimicrobial"].fillna("Unspecified")
    work["testing_lab"] = work["testing_lab"] if "testing_lab" in work.columns else pd.NA
    work["_valid"] = work["result"].astype("string").str.strip().str.upper().isin(["S", "I", "R"])
    keys = ["country", "organism", "antimicrobial", "year", "testing_lab"]
    grouped = (
        work.groupby(keys, dropna=False, sort=True)["_valid"].agg(["size", "sum"]).reset_index()
    )
    grouped.columns = [*keys, "n_records", "n_valid_result"]
    grouped["n_valid_result"] = grouped["n_valid_result"].astype(int)
    grouped["n_missing_or_invalid_result"] = grouped["n_records"] - grouped["n_valid_result"]
    grouped["pct_missing_or_invalid"] = (
        100 * grouped["n_missing_or_invalid_result"] / grouped["n_records"]
    ).round(4)
    return grouped[MISSINGNESS_COLUMNS]
