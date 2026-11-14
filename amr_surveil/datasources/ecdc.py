"""ECDC Surveillance Atlas (EARS-Net) CSV export -> aggregate surveillance context.

Counts-only, aggregate and kept as a separate stream (never joined to isolates). The importer finds
columns by name, extracts "total tested" and "resistant" counts per country/year/organism/drug group,
and **skips rows that only carry percentages** (reconstructing counts from rounded percentages would
invent precision). ``detected_columns`` is returned so the user can verify the interpretation.
"""

from __future__ import annotations

import re

import pandas as pd

from amr_surveil.constants import AGGREGATE_COLUMNS
from amr_surveil.datasources.common import PublicDataset, clean, header_key, pick
from amr_surveil.datasources.fetch import FetchError
from amr_surveil.datasources.registry import get_source

METRIC = "Resistant (R) isolates / total tested isolates, as published by ECDC (EARS-Net aggregate)"


def _parse_count(raw: object) -> float:
    """Whole-number count from text such as ``1234`` or ``1,234``; NaN for anything else."""
    text = clean(raw).replace("\u00a0", " ").strip()
    if re.fullmatch(r"\d{1,3}([ ,]\d{3})+", text):
        text = re.sub(r"[ ,]", "", text)
    if re.fullmatch(r"\d+(\.0+)?", text):
        return float(text)
    return float("nan")


def parse_ecdc_atlas(table: pd.DataFrame, *, retrieved_at: str = "") -> PublicDataset:
    keys = {header_key(c) for c in table.columns}
    needed = {"indicator", "time", "numvalue"}
    if not needed <= keys:
        raise FetchError(
            "Not an ECDC Atlas export: expected columns Indicator, Time and NumValue "
            f"(found {list(table.columns)[:10]})."
        )
    if table.empty:
        return PublicDataset(
            "ecdc_atlas",
            pd.DataFrame(),
            pd.DataFrame(),
            {"source_id": "ecdc_atlas"},
            ["The file has no data rows."],
            aggregate=pd.DataFrame(columns=AGGREGATE_COLUMNS),
        )
    population = pick(table, "Population").map(clean)
    indicator = pick(table, "Indicator").map(clean)
    unit = pick(table, "Unit").map(clean)
    country = pick(table, "RegionName", "Region name", "Country").map(clean)
    year = pd.to_numeric(pick(table, "Time").map(lambda v: clean(v)[:4]), errors="coerce")
    value = pick(table, "NumValue").map(_parse_count)

    split = population.map(lambda v: [p.strip() for p in re.split(r"\|", v, maxsplit=1)])
    organism = split.map(lambda parts: clean(parts[0]) if parts else "")
    drug = split.map(lambda parts: clean(parts[1]) if len(parts) > 1 else "")
    is_pct = indicator.str.contains(
        r"%|percent|proportion", case=False, regex=True
    ) | unit.str.contains(r"%|percent", case=False, regex=True)
    kind = pd.Series("", index=table.index)
    kind[indicator.str.contains(r"total\s+tested", case=False, regex=True) & ~is_pct] = "tested"
    kind[
        indicator.str.contains(r"^\s*R\b|resistant", case=False, regex=True)
        & ~indicator.str.contains(r"non|susceptible\b", case=False, regex=True)
        & ~is_pct
    ] = "resistant"
    work = pd.DataFrame(
        {
            "year": year,
            "country": country,
            "organism": organism,
            "antimicrobial": drug,
            "kind": kind,
            "value": value,
        }
    )
    notes: list[str] = []
    counted = work[work["kind"] != ""]
    bad_value = int(counted["value"].isna().sum())
    if bad_value:
        notes.append(
            f"{bad_value} count row(s) skipped: NumValue is not a whole number (blank, decimal or unreadable)."
        )
    work = work[
        (work["kind"] != "")
        & work["year"].notna()
        & work["value"].notna()
        & (work["country"] != "")
        & (work["antimicrobial"] != "")
    ]
    key = ["year", "country", "organism", "antimicrobial", "kind"]
    spread = work.groupby(key)["value"].transform("nunique")
    if (spread > 1).any():
        n_conflict = int(work.loc[spread > 1, key].drop_duplicates().shape[0])
        notes.append(
            f"{n_conflict} country-year-drug count(s) appear with conflicting values and were skipped "
            "(not silently picking one)."
        )
        work = work[spread <= 1]
    work = work.drop_duplicates(key)
    n_pct_only = int(is_pct.sum())
    if n_pct_only:
        notes.append(f"{n_pct_only} percentage row(s) ignored; only published counts are used.")
    wide = work.pivot_table(
        index=["year", "country", "organism", "antimicrobial"],
        columns="kind",
        values="value",
        aggfunc="first",
    ).reset_index()
    for col in ("tested", "resistant"):
        if col not in wide.columns:
            wide[col] = pd.NA
    incomplete = wide["tested"].isna() | wide["resistant"].isna()
    if incomplete.any():
        notes.append(
            f"{int(incomplete.sum())} country-year group(s) lacked either the tested or resistant count and were skipped."
        )
    wide = wide[~incomplete]
    invalid = wide["resistant"] > wide["tested"]
    if invalid.any():
        notes.append(f"{int(invalid.sum())} group(s) with resistant > tested were skipped.")
    wide = wide[~invalid & (wide["tested"] > 0)]
    info = get_source("ecdc_atlas")
    aggregate = pd.DataFrame(
        {
            "year": wide["year"].astype(int),
            "country": wide["country"],
            "source_sector": "human",
            "organism": wide["organism"],
            "antimicrobial": wide["antimicrobial"],
            "resistant_n": wide["resistant"].astype(int),
            "tested_n": wide["tested"].astype(int),
            "source_name": "ECDC Surveillance Atlas of Infectious Diseases (EARS-Net)",
            "metric_definition": METRIC,
            "source_url": info.homepage,
        }
    )[AGGREGATE_COLUMNS].reset_index(drop=True)
    if aggregate.empty:
        notes.append(
            "No usable count rows were found; check that you exported counts (not only percentages)."
        )
    return PublicDataset(
        "ecdc_atlas",
        pd.DataFrame(),
        pd.DataFrame(),
        {
            "source_id": "ecdc_atlas",
            "source_name": info.name,
            "terms": info.terms,
            "terms_url": info.terms_url,
            "terms_checked": info.terms_checked,
            "attribution": info.attribution,
            "retrieved_at": retrieved_at,
            "live_verified": False,
            "detected_columns": {"indicator": "Indicator", "year": "Time", "value": "NumValue"},
        },
        notes,
        aggregate=aggregate,
    )
