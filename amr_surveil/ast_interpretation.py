"""Optional re-check of recorded S/I/R categories against user-supplied, dated breakpoints.

No breakpoint table is embedded: EUCAST/CLSI tables are versioned, copyrighted documents that
change every year. The analyst supplies the exact table version used (``--breakpoints``); each
AST record is matched on organism, antimicrobial, standard and version. Records without a
matching breakpoint row are reported as ``no_breakpoint`` rather than guessed.

Rules (EUCAST-style, MIC in mg/L, zone in mm):
  MIC:  S if MIC <= S_le;  R if MIC > R_gt;  otherwise I (only when S_le < R_gt).
  Zone: S if zone >= S_ge; R if zone < R_lt; otherwise I (only when S_ge > R_lt).
Censored MICs are resolved only when the censoring cannot change the category; otherwise the
status is ``indeterminate``.

Units: breakpoints are in mg/L. A MIC is used only if ``mic_unit`` is blank (assumed mg/L) or an
exact equivalent (``mg/L``, ``ug/mL``, ``µg/mL``); any other unit yields ``unit_not_mg_per_l``
(or falls back to a zone diameter when one is given). Records whose standard or version is
missing never match a breakpoint row: "unknown" is not a value that two records share.
"""

from __future__ import annotations

import math

import pandas as pd

from amr_surveil.normalize import label_key, parse_mic_series

RECHECK_COLUMNS = [
    "sample_id",
    "organism",
    "antimicrobial",
    "recorded_result",
    "derived_result",
    "basis",
    "mic",
    "zone_mm",
    "breakpoint_standard",
    "breakpoint_version",
    "status",
]


_MG_PER_L = {"mg/l", "ug/ml", "µg/ml", "μg/ml", "mcg/ml"}  # numerically identical units


def _mic_unit_ok(unit: object) -> bool:
    if unit is None or pd.isna(unit):
        return True
    cleaned = str(unit).strip().lower().replace(" ", "")
    return cleaned == "" or cleaned in _MG_PER_L  # blank = assumed mg/L, as documented


def classify_mic(operator: str, value: float, s_le: float, r_gt: float) -> str | None:
    """Category from a (possibly censored) MIC, or ``None`` if indeterminate."""
    if math.isnan(value):
        return None
    if operator in ("=", ""):
        return "S" if value <= s_le else ("R" if value > r_gt else "I")
    if operator == "<=":
        return "S" if value <= s_le else None
    if operator == "<":
        return "S" if value <= s_le else None
    if operator == ">":
        return "R" if value >= r_gt else None
    if operator == ">=":
        return "R" if value > r_gt else None
    return None


def classify_zone(zone: float, s_ge: float, r_lt: float) -> str | None:
    if math.isnan(zone):
        return None
    return "S" if zone >= s_ge else ("R" if zone < r_lt else "I")


def recheck_ast_categories(ast: pd.DataFrame, breakpoints: pd.DataFrame) -> pd.DataFrame:
    """Compare recorded categories with those derived from MIC/zone and the supplied breakpoints."""
    if ast.empty or breakpoints.empty:
        return pd.DataFrame(columns=RECHECK_COLUMNS)
    bp = breakpoints.copy()
    bp["_o"] = bp["organism"].map(label_key)
    bp["_a"] = bp["antimicrobial"].map(label_key)
    bp["_s"] = bp["standard"].map(label_key)
    bp["_v"] = bp["version"].astype("string").str.strip()
    for col in ("mic_s_le", "mic_r_gt", "zone_s_ge", "zone_r_lt"):
        bp[col] = pd.to_numeric(bp.get(col), errors="coerce")
    # Missing keys must never match each other (pandas merges NaN with NaN): give breakpoint rows
    # and AST rows different sentinels.
    for col in ("_o", "_a", "_s", "_v"):
        bp[col] = bp[col].fillna("\x00missing-breakpoint-key")
    bp = bp.drop_duplicates(
        ["_o", "_a", "_s", "_v"], keep=False
    )  # ambiguous rows are never guessed

    work = ast.copy()
    for col in (
        "organism",
        "mic",
        "mic_unit",
        "zone_mm",
        "breakpoint_standard",
        "breakpoint_version",
    ):
        if col not in work.columns:
            work[col] = pd.NA
    work["_o"] = work["organism"].map(label_key)
    work["_a"] = work["antimicrobial"].map(label_key)
    work["_s"] = work["breakpoint_standard"].map(label_key)
    work["_v"] = work["breakpoint_version"].astype("string").str.strip()
    for col in ("_o", "_a", "_s", "_v"):
        work[col] = work[col].fillna("\x00missing-record-key")
    work["unit_ok"] = work["mic_unit"].map(_mic_unit_ok)
    mic = parse_mic_series(work["mic"])
    work = pd.concat([work, mic], axis=1)
    work["zone_val"] = pd.to_numeric(work["zone_mm"], errors="coerce")
    joined = work.merge(
        bp[["_o", "_a", "_s", "_v", "mic_s_le", "mic_r_gt", "zone_s_ge", "zone_r_lt"]],
        on=["_o", "_a", "_s", "_v"],
        how="left",
    )
    for col in ("mic_s_le", "mic_r_gt", "zone_s_ge", "zone_r_lt", "zone_val", "mic_value"):
        joined[col] = pd.to_numeric(joined[col], errors="coerce").astype("float64")
    joined["mic_valid"] = joined["mic_valid"].fillna(False).astype(bool)
    joined["mic_operator"] = joined["mic_operator"].fillna("").astype(str)
    derived: list[str | None] = []
    basis: list[str] = []
    status: list[str] = []
    for r in joined.itertuples(index=False):
        mic_measured = bool(r.mic_valid)
        have_mic = (
            mic_measured
            and bool(r.unit_ok)
            and not math.isnan(r.mic_s_le)
            and not math.isnan(r.mic_r_gt)
        )
        have_zone = (
            (not math.isnan(r.zone_val))
            and not math.isnan(r.zone_s_ge)
            and not math.isnan(r.zone_r_lt)
        )
        no_bp = math.isnan(r.mic_s_le) and math.isnan(r.zone_s_ge)
        if no_bp:
            derived.append(None)
            basis.append("none")
            status.append("no_breakpoint")
            continue
        if have_mic:
            cat = classify_mic(r.mic_operator, r.mic_value, r.mic_s_le, r.mic_r_gt)
            used = "mic"
        elif have_zone:
            cat = classify_zone(r.zone_val, r.zone_s_ge, r.zone_r_lt)
            used = "zone"
        else:
            derived.append(None)
            basis.append("none")
            status.append(
                "unit_not_mg_per_l" if mic_measured and not r.unit_ok else "no_measurement"
            )
            continue
        derived.append(cat)
        basis.append(used)
        recorded = str(r.result).strip().upper() if pd.notna(r.result) else ""
        if cat is None:
            status.append("indeterminate")
        else:
            status.append("concordant" if recorded == cat else "discordant")
    out = pd.DataFrame(
        {
            "sample_id": joined["sample_id"],
            "organism": joined["organism"],
            "antimicrobial": joined["antimicrobial"],
            "recorded_result": joined["result"].astype("string").str.strip().str.upper(),
            "derived_result": pd.Series(derived, dtype="string"),
            "basis": basis,
            "mic": joined["mic"],
            "zone_mm": joined["zone_mm"],
            "breakpoint_standard": joined["breakpoint_standard"],
            "breakpoint_version": joined["breakpoint_version"],
            "status": status,
        }
    )
    return out[RECHECK_COLUMNS]
