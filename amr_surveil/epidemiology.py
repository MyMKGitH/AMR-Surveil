"""Descriptive AMR indicators with explicit denominators and uncertainty."""

from __future__ import annotations

import math
from statistics import NormalDist

import pandas as pd

from amr_surveil.exceptions import ConfigurationError
from amr_surveil.exclusions import ExclusionLog
from amr_surveil.manifest import ambiguous_sample_ids, unique_manifest, unique_sample_ids
from amr_surveil.normalize import label_key, parse_iso_dates

INDICATOR_GROUP_COLUMNS = [
    "year",
    "country",
    "source_sector",
    "organism",
    "antimicrobial",
    "method",
    "breakpoint_standard",
    "breakpoint_version",
]
INDICATOR_COLUMNS = INDICATOR_GROUP_COLUMNS + [
    "n_tested",
    "n_resistant",
    "n_S",
    "n_I",
    "resistant_pct",
    "wilson_low_pct",
    "wilson_high_pct",
    "small_denominator_warning",
]
DEDUP_LOG_COLUMNS = [
    "sample_id",
    "patient_id",
    "organism",
    "specimen",
    "year",
    "group_size",
    "decision",
    "reason",
    "date_tie",
]


def wilson_interval(
    resistant: int, tested: int, confidence_level: float = 0.95
) -> tuple[float, float]:
    """Return a Wilson score interval for a binomial proportion, in percent."""
    if tested <= 0 or resistant < 0 or resistant > tested:
        raise ValueError("Require 0 <= resistant <= tested and tested > 0")
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must be between 0 and 1")
    z = NormalDist().inv_cdf(1 - (1 - confidence_level) / 2)
    p = resistant / tested
    denom = 1 + z * z / tested
    centre = (p + z * z / (2 * tested)) / denom
    half = z * math.sqrt((p * (1 - p) / tested) + (z * z / (4 * tested * tested))) / denom
    low, high = (centre - half) * 100, (centre + half) * 100
    # Exact boundary cases: avoid floating-point residue such as 2.8e-15 for 0/n.
    if resistant == 0:
        low = 0.0
    if resistant == tested:
        high = 100.0
    return max(0.0, low), min(100.0, high)


# ----------------------------------------------------------------- deduplication
def deduplicate_isolates(
    samples: pd.DataFrame,
    rule: str,
    exclusions: ExclusionLog | None = None,
    candidates: set[str] | None = None,
) -> tuple[set[str] | None, pd.DataFrame]:
    """Apply an explicitly selected first-isolate rule at *sample* level.

    ``candidates`` restricts the rule to those sample IDs (the pipeline passes the isolates that
    have at least one eligible AST record, so an isolate that was never tested cannot "win" and
    wipe out a tested repeat isolate). Configuration is still validated on the full manifest.

    Returns ``(kept_ids, log)``. ``kept_ids`` is ``None`` for rule ``none`` (keep all). For
    the other rules, every unambiguous manifest sample appears in ``log`` with a decision:
    ``kept``, ``removed`` or ``not_assignable`` (insufficient metadata; retained but not
    deduplicated). When several isolates share the earliest collection date, the lowest
    ``sample_id`` is kept and ``date_tie`` is flagged, because "first" is then arbitrary.

    The rule must follow the applicable surveillance protocol; this function does not
    choose one for you.
    """
    if rule == "none":
        return None, pd.DataFrame(columns=DEDUP_LOG_COLUMNS)
    specimen_rule = rule == "first_per_patient_species_specimen_year"
    if rule not in {"first_per_patient_species_year", "first_per_patient_species_specimen_year"}:
        raise ConfigurationError(f"Unsupported deduplication rule: {rule}")
    required = {"sample_id", "patient_id", "organism", "collection_date"} | (
        {"specimen"} if specimen_rule else set()
    )
    missing = sorted(required - set(samples.columns))
    if missing:
        raise ConfigurationError(f"{rule} requires columns in samples.csv: " + ", ".join(missing))
    manifest = unique_manifest(samples)
    if not manifest["patient_id"].notna().any():
        raise ConfigurationError(
            f"{rule} requires at least one populated patient_id; no isolates were silently deduplicated."
        )
    if candidates is not None:
        manifest = manifest[manifest["sample_id"].astype(str).isin(candidates)]
    dates = parse_iso_dates(manifest["collection_date"])
    meta = pd.DataFrame(
        {
            "sample_id": manifest["sample_id"].astype(str).to_numpy(),
            "patient_id": manifest["patient_id"].to_numpy(),
            "organism": manifest["organism"].to_numpy(),
            "specimen": manifest["specimen"].to_numpy()
            if "specimen" in manifest.columns
            else pd.NA,
            "date": dates["date"].to_numpy(),
            "precision": dates["precision"].to_numpy(),
            "year": dates["year"].to_numpy(),
        }
    )
    meta["_p"] = meta["patient_id"].map(label_key)
    meta["_o"] = meta["organism"].map(label_key)
    meta["_s"] = meta["specimen"].map(label_key)
    key_cols = ["_p", "_o", "year"] + (["_s"] if specimen_rule else [])
    eligible = meta[key_cols].notna().all(axis=1)

    rows: list[dict[str, object]] = []
    for row in meta.loc[~eligible].itertuples(index=False):
        rows.append(
            {
                "sample_id": row.sample_id,
                "patient_id": row.patient_id,
                "organism": row.organism,
                "specimen": row.specimen,
                "year": row.year,
                "group_size": 1,
                "decision": "not_assignable",
                "reason": "missing patient_id, organism, valid collection date"
                + (" or specimen" if specimen_rule else ""),
                "date_tie": False,
            }
        )
    kept: set[str] = set(meta.loc[~eligible, "sample_id"])
    ordered = meta.loc[eligible].sort_values(["date", "sample_id"], kind="stable")
    for _, group in ordered.groupby(key_cols, sort=True, dropna=False):
        first = group.iloc[0]
        # "First" is only well defined if every date in the group is a full calendar day:
        # a year- or month-only date is stored as the first day of its period and would
        # otherwise always beat a real date in that period.
        imprecise = len(group) > 1 and bool((group["precision"] != "day").any())
        tie = bool((group["date"] == first["date"]).sum() > 1) or imprecise
        kept.add(str(first["sample_id"]))
        for position, row in enumerate(group.itertuples(index=False)):
            is_first = position == 0
            rows.append(
                {
                    "sample_id": row.sample_id,
                    "patient_id": row.patient_id,
                    "organism": row.organism,
                    "specimen": row.specimen,
                    "year": row.year,
                    "group_size": len(group),
                    "decision": "kept" if is_first else "removed",
                    "reason": (
                        "first isolate in group"
                        + ("; order uncertain (imprecise collection date)" if imprecise else "")
                    )
                    if is_first
                    else f"later isolate; rule={rule}",
                    "date_tie": tie
                    if is_first or imprecise or row.date == first["date"]
                    else False,
                }
            )
            if not is_first and exclusions is not None:
                exclusions.add(
                    "deduplication",
                    "samples",
                    row.sample_id,
                    f"deduplicated:{rule}",
                    f"first isolate in group was {first['sample_id']}",
                )
    log = (
        pd.DataFrame(rows, columns=DEDUP_LOG_COLUMNS)
        .sort_values("sample_id", kind="stable")
        .reset_index(drop=True)
    )
    return kept, log


def apply_deduplication(ast: pd.DataFrame, samples: pd.DataFrame, rule: str) -> pd.DataFrame:
    """Return the AST records belonging to isolates retained by ``rule`` (compat wrapper)."""
    if rule == "none" or ast.empty:
        return ast.copy()
    kept, _ = deduplicate_isolates(samples, rule)
    assert kept is not None
    return ast[ast["sample_id"].astype(str).isin(kept)].copy()


# ------------------------------------------------------------ analysis frame
def prepare_analysis_frame(
    samples: pd.DataFrame,
    ast: pd.DataFrame,
    *,
    deduplication_rule: str = "none",
    exclusions: ExclusionLog | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Join AST to the manifest and apply every eligibility rule, logging each exclusion.

    Returns ``(frame, dedup_log)``. ``frame`` has one row per isolate/antimicrobial/method/
    breakpoint-context with a valid S/I/R category, an unambiguous sample key, a consistent
    organism and no conflicting duplicate.
    """
    log = exclusions if exclusions is not None else ExclusionLog()
    empty = pd.DataFrame(
        columns=["sample_id", *INDICATOR_GROUP_COLUMNS, "result", "testing_lab", "specimen"]
    )
    if ast.empty or "sample_id" not in samples.columns or samples.empty:
        return empty, pd.DataFrame(columns=DEDUP_LOG_COLUMNS)

    work = ast.copy()
    if "antimicrobial" not in work.columns:
        work["antimicrobial"] = pd.NA
    work["result"] = work["result"].astype("string").str.strip().str.upper()
    record_id = (
        work["sample_id"].astype("string").fillna("<missing>")
        + "|"
        + work["antimicrobial"].astype("string").fillna("<missing>")
    )
    bad_result = ~work["result"].isin(["S", "I", "R"])
    for rid in record_id[bad_result]:
        log.add("indicators", "ast", rid, "missing_or_invalid_result", "no valid S/I/R category")
    work = work[~bad_result].copy()
    no_drug = work["antimicrobial"].isna()
    for rid in work.loc[no_drug, "sample_id"].astype("string").fillna("<missing>") + "|<missing>":
        log.add("indicators", "ast", rid, "missing_antimicrobial", "record has no antimicrobial")
    work = work[~no_drug].copy()

    ambiguous = ambiguous_sample_ids(samples)
    unique_ids = unique_sample_ids(samples)
    sid = work["sample_id"].astype("string")
    for rid in sid[sid.isin(ambiguous)].dropna().unique():
        log.add(
            "indicators",
            "ast",
            rid,
            "ambiguous_sample_id",
            "sample_id is duplicated in the manifest",
        )
    unknown = ~sid.isin(unique_ids) & ~sid.isin(ambiguous)
    for rid in sid[unknown].dropna().unique():
        log.add(
            "indicators", "ast", rid, "unmatched_sample_id", "sample_id absent from the manifest"
        )
    work = work[sid.isin(unique_ids)].copy()

    clean_samples = unique_manifest(samples)

    meta_cols = [
        c
        for c in [
            "sample_id",
            "organism",
            "collection_date",
            "country",
            "source_sector",
            "specimen",
        ]
        if c in clean_samples.columns
    ]
    merged = work.merge(
        clean_samples[meta_cols],
        on="sample_id",
        how="inner",
        suffixes=("_ast", "_sample"),
        validate="many_to_one",
    )
    if "organism_ast" in merged.columns:
        ast_org = merged["organism_ast"].map(label_key)
        sample_org = merged["organism_sample"].map(label_key)
        mismatch = ast_org.notna() & sample_org.notna() & (ast_org != sample_org)
        for rid in merged.loc[mismatch, "sample_id"].unique():
            log.add(
                "indicators",
                "ast",
                rid,
                "organism_mismatch",
                "AST organism conflicts with the manifest",
            )
        merged = merged[~mismatch].copy()
        merged["organism"] = merged["organism_ast"].fillna(merged["organism_sample"])
    elif "organism" not in merged.columns:
        merged["organism"] = pd.NA
    # Deduplicate only among isolates that actually contribute an eligible AST record, so an
    # untested (or wholly invalid) first isolate cannot remove a tested repeat isolate.
    kept, dedup_log = deduplicate_isolates(
        clean_samples,
        deduplication_rule,
        log,
        candidates=set(merged["sample_id"].astype(str)),
    )
    if kept is not None:
        merged = merged[merged["sample_id"].astype(str).isin(kept)].copy()
    for col in ("country", "source_sector"):
        if col not in merged.columns:
            merged[col] = "Unknown"
        merged[col] = merged[col].fillna("Unknown")
    if "specimen" not in merged.columns:
        merged["specimen"] = pd.NA
    if "collection_date" in merged.columns:
        merged["year"] = parse_iso_dates(merged["collection_date"])["year"].astype("string")
    else:
        merged["year"] = pd.Series(pd.NA, index=merged.index, dtype="string")
    for col in ("method", "breakpoint_standard", "breakpoint_version"):
        if col not in merged.columns:
            merged[col] = "Unspecified"
        merged[col] = merged[col].fillna("Unspecified")
    if "testing_lab" not in merged.columns:
        merged["testing_lab"] = pd.NA

    key = [
        "sample_id",
        "year",
        "country",
        "source_sector",
        "organism",
        "antimicrobial",
        "method",
        "breakpoint_standard",
        "breakpoint_version",
    ]
    before = len(merged)
    merged = merged.drop_duplicates(subset=key + ["result"])
    n_collapsed = before - len(merged)
    if n_collapsed:
        log.add(
            "indicators",
            "ast",
            f"{n_collapsed} row(s)",
            "exact_duplicate_collapsed",
            "identical isolate/drug/method/breakpoint/category rows counted once",
        )
    conflict_counts = merged.groupby(key, dropna=False)["result"].transform("nunique")
    conflicting = conflict_counts > 1
    for rid in (
        merged.loc[conflicting, "sample_id"].astype(str)
        + "|"
        + merged.loc[conflicting, "antimicrobial"].astype(str)
    ).unique():
        log.add(
            "indicators",
            "ast",
            rid,
            "conflicting_duplicate_ast",
            "conflicting S/I/R for the same isolate/drug/method/breakpoint",
        )
    merged = merged[~conflicting].copy()
    keep_cols = ["sample_id", *INDICATOR_GROUP_COLUMNS, "result", "testing_lab", "specimen"]
    keep_cols = list(dict.fromkeys(keep_cols))
    return merged[keep_cols].reset_index(drop=True), dedup_log


def indicators_from_frame(
    frame: pd.DataFrame,
    *,
    confidence_level: float = 0.95,
    low_denominator_threshold: int = 10,
) -> pd.DataFrame:
    """R/(S+I+R) by group from a prepared analysis frame; percentages are never averaged."""
    if frame.empty:
        return pd.DataFrame(columns=INDICATOR_COLUMNS)
    records: list[dict[str, object]] = []
    for group_key, group in frame.groupby(INDICATOR_GROUP_COLUMNS, dropna=False, sort=True):
        per_isolate = group.drop_duplicates("sample_id")
        tested = int(per_isolate["sample_id"].nunique())
        if tested == 0:
            continue
        counts = per_isolate["result"].value_counts()
        resistant = int(counts.get("R", 0))
        low, high = wilson_interval(resistant, tested, confidence_level)
        record = dict(zip(INDICATOR_GROUP_COLUMNS, group_key, strict=True))
        record.update(
            {
                "n_tested": tested,
                "n_resistant": resistant,
                "n_S": int(counts.get("S", 0)),
                "n_I": int(counts.get("I", 0)),
                "resistant_pct": round(100 * resistant / tested, 4),
                "wilson_low_pct": round(low, 4),
                "wilson_high_pct": round(high, 4),
                "small_denominator_warning": tested < low_denominator_threshold,
            }
        )
        records.append(record)
    return pd.DataFrame(records, columns=INDICATOR_COLUMNS)


def resistance_indicators(
    samples: pd.DataFrame,
    ast: pd.DataFrame,
    *,
    deduplication_rule: str = "none",
    confidence_level: float = 0.95,
    low_denominator_threshold: int = 10,
    exclusions: ExclusionLog | None = None,
) -> pd.DataFrame:
    """Compute R/(S+I+R) by year, sector and method/breakpoint context (public convenience API)."""
    frame, _ = prepare_analysis_frame(
        samples, ast, deduplication_rule=deduplication_rule, exclusions=exclusions
    )
    return indicators_from_frame(
        frame,
        confidence_level=confidence_level,
        low_denominator_threshold=low_denominator_threshold,
    )


def aggregate_surveillance_context(aggregate: pd.DataFrame) -> pd.DataFrame:
    """Calculate percentages from supplied counts while retaining source rows separately."""
    if aggregate.empty:
        return pd.DataFrame(columns=[*aggregate.columns, "resistant_pct", "evidence_level"])
    frame = aggregate.copy()
    frame["resistant_n"] = pd.to_numeric(frame["resistant_n"], errors="coerce")
    frame["tested_n"] = pd.to_numeric(frame["tested_n"], errors="coerce")
    valid = (
        frame["resistant_n"].notna()
        & frame["tested_n"].notna()
        & (frame["tested_n"] > 0)
        & (frame["resistant_n"] >= 0)
        & (frame["resistant_n"] <= frame["tested_n"])
    )
    frame = frame[valid].copy()
    frame["resistant_pct"] = (100 * frame["resistant_n"] / frame["tested_n"]).round(4)
    frame["evidence_level"] = "aggregate surveillance; not isolate-linked"
    return frame
