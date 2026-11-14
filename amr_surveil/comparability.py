"""Threats to comparability across years, laboratories and methods.

For each organism/antimicrobial/country/sector series the module raises *flags*; it never
adjusts, smooths or models the data. A flag means "do not read a trend into these numbers until
the cause is understood".
"""

from __future__ import annotations

import pandas as pd

SERIES = ["country", "source_sector", "organism", "antimicrobial"]
FLAG_COLUMNS = [*SERIES, "flag", "severity", "detail"]


def _fmt(items) -> str:
    return ", ".join(sorted(map(str, items)))


def _emit_if_varying(emit, key, known: pd.DataFrame, col: str, flag: str, label: str) -> None:
    """Flag a series whose set of values for ``col`` differs between years."""
    sets = known.groupby("_y")[col].apply(lambda s: frozenset(s.dropna().astype(str)))
    sets = sets[sets.map(len) > 0]
    if len(sets) >= 2 and len(set(sets)) > 1:
        emit(
            key,
            flag,
            "WARNING",
            f"{label} differs between years: "
            + "; ".join(f"{y}: {_fmt(v)}" for y, v in sets.items()),
        )


def comparability_flags(
    frame: pd.DataFrame,
    *,
    small_n: int = 10,
    volume_jump_ratio: float = 2.0,
) -> pd.DataFrame:
    """Return one row per (series, flag). ``frame`` is the prepared per-isolate analysis frame."""
    if frame.empty:
        return pd.DataFrame(columns=FLAG_COLUMNS)
    out: list[dict[str, object]] = []

    def emit(key, flag, severity, detail):
        out.append(
            {
                **dict(zip(SERIES, key, strict=True)),
                "flag": flag,
                "severity": severity,
                "detail": detail,
            }
        )

    for key, g in frame.groupby(SERIES, dropna=False, sort=True):
        n_unknown = int(g["year"].isna().sum())
        if n_unknown:
            emit(
                key,
                "unknown_year_records",
                "WARNING",
                f"{n_unknown} record(s) lack a valid collection year and cannot be placed in a trend",
            )
        known = g[g["year"].notna()].copy()
        known["_y"] = known["year"].astype(int)
        distinct_years = sorted(known["_y"].unique())
        if len(distinct_years) < 2:
            emit(
                key,
                "single_year_only",
                "INFO",
                "fewer than two years: no temporal comparison is possible",
            )
            continue
        per_year_n = known.drop_duplicates(["sample_id", "_y"]).groupby("_y")["sample_id"].nunique()

        small = per_year_n[per_year_n < small_n]
        if not small.empty:
            emit(
                key,
                "small_denominator_years",
                "WARNING",
                f"years with fewer than {small_n} isolates: {_fmt(f'{y} (n={n})' for y, n in small.items())}",
            )
        span = set(range(min(distinct_years), max(distinct_years) + 1))
        gaps = span - set(distinct_years)
        if gaps:
            emit(key, "missing_years", "WARNING", f"no isolates in year(s): {_fmt(gaps)}")
        for prev, cur in zip(distinct_years, distinct_years[1:], strict=False):
            if cur - prev != 1:
                continue
            a, b = int(per_year_n[prev]), int(per_year_n[cur])
            if a and b and (b / a >= volume_jump_ratio or a / b >= volume_jump_ratio):
                emit(
                    key,
                    "testing_volume_jump",
                    "WARNING",
                    f"tested isolates changed from {a} ({prev}) to {b} ({cur}); possible change in panel, "
                    "laboratory participation or sampling",
                )

        known["_bp"] = (
            known["breakpoint_standard"].astype(str) + " " + known["breakpoint_version"].astype(str)
        )
        for column, flag, label in (
            ("_bp", "breakpoint_version_changes_across_years", "breakpoint standard/version"),
            ("method", "ast_method_changes_across_years", "AST method"),
            ("testing_lab", "testing_laboratory_set_changes", "testing laboratory set"),
            ("specimen", "specimen_mix_changes", "specimen types"),
        ):
            _emit_if_varying(emit, key, known, column, flag, label)
    return pd.DataFrame(out, columns=FLAG_COLUMNS)
