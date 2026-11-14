"""NCBI AMRFinderPlus report adapter.

AMRFinderPlus is the native AMR engine inside FishAMR-Link (its Dockerfile pins v3.12.8 and
its code reads the ``Gene symbol``, ``Sequence name``, ``Class``, ``Subclass`` and ``Method``
columns). This adapter converts the tab-separated report(s) into AMR-Surveil's canonical
``amr_calls`` and ``screening`` tables.

Why a screening table matters
-----------------------------
AMRFinderPlus prints *positive* hits only. A genome with no hits yields a header-only file.
Recording each report as a screened sample is what later allows a "no mapped gene" result to be
treated as a negative rather than as missing data. Tool and database versions are not written
into the report, so supply them explicitly; they are recorded, never invented.

Both the v3.x header style (``Gene symbol``, ``Element type``) and the v4.x style
(``Element symbol``, ``Type``) are accepted.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from amr_surveil.constants import AMR_CALL_COLUMNS, SCREENING_COLUMNS
from amr_surveil.exceptions import ConfigurationError, InputFileError
from amr_surveil.io import read_csv, write_csv

ALIASES: dict[str, list[str]] = {
    "name": ["Name"],
    "contig": ["Contig id"],
    "start": ["Start"],
    "end": ["Stop", "End"],
    "gene": ["Gene symbol", "Element symbol"],
    "sequence_name": ["Sequence name", "Element name"],
    "element_type": ["Element type", "Type"],
    "klass": ["Class"],
    "subclass": ["Subclass"],
    "method": ["Method"],
    "coverage": ["% Coverage of reference sequence", "% Coverage of reference"],
    "identity": ["% Identity to reference sequence", "% Identity to reference"],
}
_BLANK = {"", "NA", "N/A", "-", "nan"}


def _pick(columns: list[str], key: str) -> str | None:
    for candidate in ALIASES[key]:
        if candidate in columns:
            return candidate
    return None


def _col_or_na(raw: pd.DataFrame, cols: list[str], key: str) -> pd.Series:
    col = _pick(cols, key)
    if col is None:
        return pd.Series(pd.NA, index=raw.index, dtype="string")
    return _clean(raw[col])


def _clean(series: pd.Series) -> pd.Series:
    s = series.astype("string").str.strip()
    return s.mask(s.isin(_BLANK), pd.NA)


def normalize_amrfinder_reports(
    inputs: list[str | Path],
    output_calls: str | Path,
    output_screening: str | Path | None = None,
    *,
    sample_ids: list[str] | None = None,
    organism: str | None = None,
    tool_version: str | None = None,
    database_version: str | None = None,
    database_name: str = "NCBI AMRFinderPlus database",
    element_types: tuple[str, ...] = ("AMR",),
    identity_threshold_pct: float | None = None,
    coverage_threshold_pct: float | None = None,
) -> dict[str, object]:
    """Normalise one or more AMRFinderPlus TSV reports.

    Sample IDs come from the report's ``Name`` column when present, otherwise from ``sample_ids``
    (same order as ``inputs``), otherwise from the file name stem. Returns a summary dictionary.
    """
    if not inputs:
        raise ConfigurationError("At least one AMRFinderPlus report is required")
    if sample_ids is not None and len(sample_ids) != len(inputs):
        raise ConfigurationError("sample_ids must have the same length as inputs")
    wanted = {e.upper() for e in element_types}
    call_frames: list[pd.DataFrame] = []
    screened: dict[str, str] = {}
    excluded: dict[str, int] = {}
    n_rows = 0
    for position, path in enumerate(inputs):
        path = Path(path)
        if path.exists() and path.stat().st_size == 0:
            raise InputFileError(
                f"{path}: file is empty (a report with no hits still has a header row)"
            )
        raw = read_csv(path, table=f"amrfinder:{path.name}", sep="\t")
        cols = list(raw.columns)
        gene_col = _pick(cols, "gene")
        if gene_col is None:
            raise InputFileError(
                f"{path}: no gene column found (expected one of {ALIASES['gene']}); is this an AMRFinderPlus report?"
            )
        fallback_id = sample_ids[position] if sample_ids is not None else path.stem
        name_col = _pick(cols, "name")
        if name_col is not None and raw[name_col].notna().any():
            sid = raw[name_col].fillna(fallback_id).astype(str)
        else:
            sid = pd.Series(fallback_id, index=raw.index, dtype="string")
        file_samples = set(sid) if len(raw) else {fallback_id}
        clash = file_samples & set(screened)
        if clash:
            raise ConfigurationError(
                f"Sample(s) appear in more than one report: {', '.join(sorted(clash))}"
            )
        for s in file_samples:
            screened[s] = path.name
        n_rows += len(raw)
        if raw.empty:
            continue
        gene = _clean(raw[gene_col])
        seq_name_col = _pick(cols, "sequence_name")
        if seq_name_col:
            gene = gene.fillna(_clean(raw[seq_name_col]))
        etype_col = _pick(cols, "element_type")
        etype = (
            _clean(raw[etype_col])
            if etype_col
            else pd.Series("AMR", index=raw.index, dtype="string")
        )
        keep = etype.str.upper().isin(wanted) & gene.notna()
        for label, count in etype[~keep].fillna("missing_gene").value_counts().items():
            excluded[str(label)] = excluded.get(str(label), 0) + int(count)
        subclass_col, class_col = _pick(cols, "subclass"), _pick(cols, "klass")
        det = (
            _clean(raw[subclass_col])
            if subclass_col
            else pd.Series(pd.NA, index=raw.index, dtype="string")
        )
        if class_col:
            det = det.fillna(_clean(raw[class_col]))
        frame = pd.DataFrame(
            {
                "sample_id": sid,
                "organism": organism,
                "gene": gene,
                "identity_pct": _col_or_na(raw, cols, "identity"),
                "coverage_pct": _col_or_na(raw, cols, "coverage"),
                "contig": _col_or_na(raw, cols, "contig"),
                "start": _col_or_na(raw, cols, "start"),
                "end": _col_or_na(raw, cols, "end"),
                "method": "AMRFinderPlus",
                "database": database_name,
                "database_version": database_version,
                "source_system": "AMRFinderPlus",
                "element_type": etype,
                "determinant_class": det,
                "hit_method": _col_or_na(raw, cols, "method"),
                "evidence_tier": "curated_database",
            }
        )[keep.to_numpy()]
        call_frames.append(frame)
    calls = (
        pd.concat(call_frames, ignore_index=True)
        if call_frames
        else pd.DataFrame(columns=AMR_CALL_COLUMNS)
    )
    for col in AMR_CALL_COLUMNS:
        if col not in calls.columns:
            calls[col] = pd.NA
    calls = calls[AMR_CALL_COLUMNS]
    write_csv(calls, output_calls)
    screening = pd.DataFrame(
        {
            "sample_id": sorted(screened),
            "method": "AMRFinderPlus",
            "tool_version": tool_version,
            "database": database_name,
            "database_version": database_version,
            "identity_threshold_pct": identity_threshold_pct,
            "coverage_threshold_pct": coverage_threshold_pct,
            "source_system": "AMRFinderPlus",
        }
    )[SCREENING_COLUMNS]
    if output_screening is not None:
        write_csv(screening, output_screening)
    return {
        "n_reports": len(inputs),
        "n_samples_screened": len(screened),
        "n_report_rows": n_rows,
        "n_calls_written": len(calls),
        "excluded_by_element_type": excluded,
        "tool_version_recorded": tool_version is not None,
        "database_version_recorded": database_version is not None,
    }
