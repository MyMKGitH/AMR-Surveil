"""Strict, deterministic parsing and label-harmonisation helpers.

Nothing here guesses. Dates must be ISO 8601 (``YYYY``, ``YYYY-MM`` or ``YYYY-MM-DD``);
ambiguous day/month orderings are rejected rather than inferred, because pandas format
inference makes the outcome depend on the first row of a file. MIC strings keep their
censoring operator. Label harmonisation only folds case and whitespace and always reports
the variants it merged.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

import pandas as pd

_ISO_DAY = re.compile(r"\d{4}-\d{2}-\d{2}")
_ISO_MONTH = re.compile(r"\d{4}-\d{2}")
_ISO_YEAR = re.compile(r"\d{4}")
_MIC_PATTERN = re.compile(
    r"^(<=|>=|<|>|=|≤|≥)?\s*((?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)$"
)
_OPERATOR_CANON = {"≤": "<=", "≥": ">=", "": "="}


# --------------------------------------------------------------------------- dates
def parse_iso_dates(values: pd.Series) -> pd.DataFrame:
    """Parse ISO-8601 partial dates without format inference.

    Returns a frame with ``date`` (datetime64, first day of the period), ``precision``
    (``day``/``month``/``year``), ``year`` (nullable Int64), ``non_iso`` (bool: a value was
    present but not in an accepted ISO form) and ``impossible`` (bool: ISO-shaped but not a
    real calendar date, e.g. 2024-02-31).
    """
    text = values.astype("string").str.strip()
    present = text.notna()
    is_day = text.str.fullmatch(_ISO_DAY).fillna(False).astype(bool)
    is_month = text.str.fullmatch(_ISO_MONTH).fillna(False).astype(bool)
    is_year = text.str.fullmatch(_ISO_YEAR).fillna(False).astype(bool)
    shaped = is_day | is_month | is_year

    parsed = pd.Series(pd.NaT, index=values.index, dtype="datetime64[ns]")
    precision = pd.Series(pd.NA, index=values.index, dtype="string")
    for mask, fmt, label in (
        (is_day, "%Y-%m-%d", "day"),
        (is_month, "%Y-%m", "month"),
        (is_year, "%Y", "year"),
    ):
        if mask.any():
            converted = pd.to_datetime(text[mask], format=fmt, errors="coerce")
            # datetime64[ns] holds 1677-2262 only; anything outside is not a usable date here.
            converted = converted.where(converted.dt.year.between(1678, 2261))
            parsed.loc[mask] = converted.astype("datetime64[ns]")
            precision.loc[mask] = label
    ok = parsed.notna()
    precision = precision.where(ok, pd.NA)
    year = parsed.dt.year.astype("Int64")
    return pd.DataFrame(
        {
            "date": parsed,
            "precision": precision,
            "year": year,
            "non_iso": present & ~shaped,
            "impossible": shaped & ~ok,
        }
    )


def collection_year(values: pd.Series) -> pd.Series:
    """Calendar year (nullable Int64) from strict ISO dates; NA when unparsable."""
    return parse_iso_dates(values)["year"]


# ----------------------------------------------------------------------------- MIC
@dataclass(frozen=True)
class MicValue:
    operator: str
    value: float
    valid: bool


def parse_mic(raw: object) -> MicValue:
    """Parse a MIC such as ``0.5``, ``<=0.25`` or ``>32`` (``≤``/``≥`` accepted).

    A decimal comma is deliberately *not* accepted: ``1,000`` is ambiguous. Such values are
    reported as invalid for the analyst to correct at source.
    """
    if raw is None or (isinstance(raw, float) and math.isnan(raw)) or pd.isna(raw):
        return MicValue("", math.nan, False)
    text = str(raw).strip()
    match = _MIC_PATTERN.match(text)
    if not match:
        return MicValue("", math.nan, False)
    operator = _OPERATOR_CANON.get(match.group(1) or "", match.group(1) or "=")
    return MicValue(operator, float(match.group(2)), True)


def parse_mic_series(values: pd.Series) -> pd.DataFrame:
    """Vectorised wrapper returning ``mic_operator``, ``mic_value`` and ``mic_valid``."""
    parsed = [parse_mic(v) for v in values.astype(object)]
    return pd.DataFrame(
        {
            "mic_operator": [p.operator for p in parsed],
            "mic_value": [p.value for p in parsed],
            "mic_valid": [p.valid for p in parsed],
        },
        index=values.index,
    )


# ----------------------------------------------------------------------- label keys
def label_key(value: object) -> str | None:
    """Case-/whitespace-folded comparison key (``None`` for missing)."""
    if value is None or pd.isna(value):
        return None
    return re.sub(r"\s+", " ", str(value)).strip().casefold()


def canonical_spellings(values: pd.Series) -> dict[str, str]:
    """Map each case-/whitespace-folded key to its canonical spelling.

    The canonical spelling is the most frequent variant (ties broken alphabetically), so the
    choice is deterministic. Pass values pooled from several tables to harmonise them jointly.
    """
    text = values.astype("string").str.replace(r"\s+", " ", regex=True).str.strip()
    counts = text.value_counts(dropna=True)
    canonical: dict[str, str] = {}
    for spelling, _count in sorted(counts.items(), key=lambda kv: (-kv[1], str(kv[0]))):
        key = label_key(spelling)
        if key is not None:
            canonical.setdefault(key, str(spelling))
    return canonical


def harmonise_labels(
    values: pd.Series, canonical: dict[str, str] | None = None
) -> tuple[pd.Series, list[dict[str, object]]]:
    """Fold case/whitespace variants of the same label onto one canonical spelling.

    ``canonical`` may be supplied (see :func:`canonical_spellings`) so several tables share one
    spelling; otherwise it is derived from ``values``. Returns the harmonised series and one
    record per key whose spelling was changed or that had several variants, so the caller can
    raise an audit finding. Synonyms (``E. coli`` vs ``Escherichia coli``) are *not* merged: that
    requires a curated dictionary upstream.
    """
    text = values.astype("string").str.replace(r"\s+", " ", regex=True).str.strip()
    keys = text.map(label_key)
    canonical = canonical if canonical is not None else canonical_spellings(values)
    mapped = keys.map(canonical).astype("string")
    merged = []
    for key in sorted(k for k in keys.dropna().unique()):
        seen = sorted(set(text[keys == key].dropna().astype(str)))
        if len(seen) > 1 or seen[0] != canonical[key]:
            merged.append({"canonical": canonical[key], "variants": seen})
    # Missing values stay missing.
    return mapped.where(text.notna(), pd.NA), merged
