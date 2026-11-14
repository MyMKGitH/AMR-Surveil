"""Cumulative antibiogram (percentage susceptible) in the style of CLSI M39.

Built from the pipeline's ``resistance_by_group`` table, so it inherits the configured
deduplication rule and the AST eligibility logging. Differences from the resistance proportion used
elsewhere in AMR-Surveil are deliberate and stated: an antibiogram reports **% susceptible = S /
(S+I+R)** (I counts as not-susceptible), whereas the resistance indicator is R / (S+I+R).
Sectors (human / animal / environment ...) get separate rows. Countries and years inside the chosen
filters are summed; note that with a *per-year* deduplication rule a patient can be counted once in
each year of a multi-year antibiogram (CLSI M39 asks for once per analysis period), so prefer the
``none``/protocol rule you actually need. Rows below ``min_n`` isolates are flagged (CLSI recommends at least 30 isolates per
organism-drug combination for a reliable estimate). Different contexts (method, breakpoint standard
and version) are never merged.
"""

from __future__ import annotations

import pandas as pd

from amr_surveil.epidemiology import wilson_interval

KEYS = [
    "organism",
    "antimicrobial",
    "source_sector",
    "method",
    "breakpoint_standard",
    "breakpoint_version",
]
COLUMNS = [
    *KEYS,
    "n_tested",
    "n_S",
    "n_I",
    "n_R",
    "susceptible_pct",
    "wilson_low_pct",
    "wilson_high_pct",
    "meets_min_n",
]


def cumulative_antibiogram(
    indicators: pd.DataFrame,
    *,
    year_from: int | None = None,
    year_to: int | None = None,
    source_sector: str | None = None,
    country: str | None = None,
    min_n: int = 30,
) -> pd.DataFrame:
    if min_n < 1:
        raise ValueError("min_n must be at least 1")
    if indicators.empty:
        return pd.DataFrame(columns=COLUMNS)
    work = indicators.copy()
    year = pd.to_numeric(work["year"], errors="coerce")
    if year_from is not None:
        work = work[year >= year_from]
        year = year.loc[work.index]
    if year_to is not None:
        work = work[year <= year_to]
    if source_sector is not None:
        work = work[work["source_sector"] == source_sector]
    if country is not None:
        work = work[work["country"] == country]
    if work.empty:
        return pd.DataFrame(columns=COLUMNS)
    for col in KEYS:
        work[col] = work[col].astype("string").fillna("")
    grouped = (
        work.groupby(KEYS, sort=True)[["n_tested", "n_S", "n_I", "n_resistant"]].sum().reset_index()
    )
    grouped = grouped.rename(columns={"n_resistant": "n_R"})
    rows = []
    for rec in grouped.to_dict("records"):
        tested, sus = int(rec["n_tested"]), int(rec["n_S"])
        low, high = wilson_interval(sus, tested)
        rows.append(
            {
                **{k: rec[k] for k in KEYS},
                "n_tested": tested,
                "n_S": sus,
                "n_I": int(rec["n_I"]),
                "n_R": int(rec["n_R"]),
                "susceptible_pct": round(100 * sus / tested, 2),
                "wilson_low_pct": round(low, 2),
                "wilson_high_pct": round(high, 2),
                "meets_min_n": tested >= min_n,
            }
        )
    return pd.DataFrame(rows, columns=COLUMNS)
