"""Country-level resistance tables for maps (pure pandas; the figure itself is drawn in the UI).

Guard rails, all enforced here rather than left to the user:

* **One context per map.** Counts are summed (never percentages averaged) and only inside a single
  method + breakpoint standard + version; a selection spanning several contexts raises.
* **One sector per map.** Human, animal and environmental isolates are not pooled.
* **Small denominators are suppressed**, not shaded: a country below ``min_n`` tested isolates is
  returned with ``suppressed=True`` and no colour value.
* **Unmappable country names are reported**, never dropped silently.
"""

from __future__ import annotations

import pandas as pd

from amr_surveil.datasources.geo import country_to_iso3, iso3_to_name
from amr_surveil.epidemiology import wilson_interval

MAP_COLUMNS = [
    "country",
    "iso3",
    "n_resistant",
    "n_tested",
    "resistant_pct",
    "wilson_low_pct",
    "wilson_high_pct",
    "n_years",
    "suppressed",
]
CONTEXT_COLUMNS = ["method", "breakpoint_standard", "breakpoint_version"]


def country_resistance(
    frame: pd.DataFrame,
    *,
    organism: str,
    antimicrobial: str,
    year_from: int | None = None,
    year_to: int | None = None,
    source_sector: str | None = None,
    min_n: int = 20,
    resistant_col: str = "n_resistant",
    tested_col: str = "n_tested",
) -> tuple[pd.DataFrame, list[str]]:
    """Per-country R/(S+I+R) for one organism/antimicrobial. Returns ``(table, unmapped_countries)``."""
    if min_n < 1:
        raise ValueError("min_n must be at least 1")
    if frame.empty:
        return pd.DataFrame(columns=MAP_COLUMNS), []
    work = frame.copy()
    work = work[(work["organism"] == organism) & (work["antimicrobial"] == antimicrobial)]
    year = pd.to_numeric(work["year"], errors="coerce")
    if year_from is not None:
        work = work[year >= year_from]
        year = year.loc[work.index]
    if year_to is not None:
        work = work[year <= year_to]
    if source_sector is not None and "source_sector" in work.columns:
        work = work[work["source_sector"] == source_sector]
    if work.empty:
        return pd.DataFrame(columns=MAP_COLUMNS), []
    if "source_sector" in work.columns and work["source_sector"].nunique() > 1:
        raise ValueError(
            "The selection mixes source sectors "
            f"{sorted(work['source_sector'].astype(str).unique())}; pick one sector."
        )
    present = [c for c in CONTEXT_COLUMNS if c in work.columns]
    if present:
        contexts = work[present].fillna("").astype(str).drop_duplicates()
        if len(contexts) > 1:
            raise ValueError(
                "The selection spans several AST method/breakpoint contexts "
                f"({len(contexts)}); choose one so counts are not pooled across incomparable tests."
            )
    work["country"] = work["country"].astype("string").fillna("").str.strip()
    blank = work["country"] == ""
    work = work[~blank]
    records = []
    unmapped: set[str] = set()
    # Spellings of one country ("USA", "United States") are merged by ISO code before any
    # suppression decision; unrecognised names stay separate and are reported.
    # Built from plain Python objects: ``Series.map`` / ``.where`` on string dtypes behave
    # differently across pandas versions when every code is missing.
    countries = work["country"].astype(str).tolist()
    iso_codes = [country_to_iso3(name) for name in countries]
    work["_iso3"] = pd.Series(iso_codes, index=work.index, dtype="object")
    work["_key"] = pd.Series(
        [iso if iso else f"name:{name}" for iso, name in zip(iso_codes, countries, strict=True)],
        index=work.index,
        dtype="object",
    )
    for _, group in work.groupby("_key", sort=True):
        iso3 = group["_iso3"].iloc[0]
        iso3 = None if pd.isna(iso3) else str(iso3)
        names = sorted(group["country"].unique())
        country = (iso3_to_name(iso3) if iso3 else None) or names[0]
        if iso3 is None:
            unmapped.add(str(country))
        tested = int(pd.to_numeric(group[tested_col]).sum())
        resistant = int(pd.to_numeric(group[resistant_col]).sum())
        if tested <= 0:
            continue
        low, high = wilson_interval(resistant, tested)
        suppressed = tested < min_n
        records.append(
            {
                "country": country,
                "iso3": iso3,
                "n_resistant": resistant,
                "n_tested": tested,
                "resistant_pct": None if suppressed else round(100 * resistant / tested, 2),
                "wilson_low_pct": None if suppressed else round(low, 2),
                "wilson_high_pct": None if suppressed else round(high, 2),
                "n_years": int(pd.to_numeric(group["year"], errors="coerce").nunique()),
                "suppressed": suppressed,
            }
        )
    return pd.DataFrame(records, columns=MAP_COLUMNS), sorted(unmapped)
