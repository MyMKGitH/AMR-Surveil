"""Input/output helpers with stable string handling and provenance utilities."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

from amr_surveil.exceptions import InputFileError


def read_csv(
    path: str | Path | None,
    *,
    required: list[str] | None = None,
    table: str = "input",
    sep: str = ",",
) -> pd.DataFrame:
    """Read a delimited table without silently converting identifiers or blank values.

    ``sep`` is ``","`` for CSV or ``"\\t"`` for TSV (e.g. AMRFinderPlus output).
    """
    if path is None:
        return pd.DataFrame(columns=required or [])
    csv_path = Path(path)
    if not csv_path.exists() or not csv_path.is_file():
        raise InputFileError(f"{table}: file does not exist or is not a regular file: {csv_path}")
    try:
        with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter=sep)
            try:
                raw_header = next(reader)
            except StopIteration as exc:
                raise InputFileError(f"{table}: CSV is empty and has no header row") from exc
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise InputFileError(f"{table}: could not read CSV header {csv_path}: {exc}") from exc

    header = [str(col).strip() for col in raw_header]
    if not header or any(not col for col in header):
        raise InputFileError(f"{table}: blank column names are not allowed")
    duplicates = sorted(col for col, count in Counter(header).items() if count > 1)
    if duplicates:
        raise InputFileError(
            f"{table}: duplicate column name(s) are not allowed: {', '.join(duplicates)}"
        )

    try:
        frame = pd.read_csv(
            csv_path,
            dtype=str,
            keep_default_na=False,
            na_filter=False,
            encoding="utf-8-sig",
            sep=sep,
        )
    except (OSError, UnicodeDecodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise InputFileError(f"{table}: could not read CSV {csv_path}: {exc}") from exc
    frame.columns = [str(col).strip() for col in frame.columns]
    for col in frame.columns:
        frame[col] = frame[col].astype("string").str.strip()
        frame[col] = frame[col].mask(frame[col] == "", pd.NA)
    missing = sorted(set(required or []) - set(frame.columns))
    if missing:
        raise InputFileError(f"{table}: missing required column(s): {', '.join(missing)}")
    return frame


_FORMULA_START = re.compile(r"^(?:[=+@\t\r]|-(?![0-9.]))")


def neutralise_formulas(frame: pd.DataFrame) -> pd.DataFrame:
    """Prefix text cells that a spreadsheet would execute as a formula with a single quote.

    Echoed user-supplied text (sample IDs, organism names, free-text notes) is untrusted: a cell
    starting with ``=``, ``+``, ``@``, a tab/CR, or ``-`` followed by a non-digit can run as a
    formula when the CSV is opened in Excel/LibreOffice (CSV injection, OWASP). Numbers and
    negative numeric text such as ``-1.5`` are left untouched.
    """
    out = frame.copy()
    for pos in range(out.shape[1]):  # positional: tolerates duplicate column labels
        series = out.iloc[:, pos]
        if not (
            pd.api.types.is_object_dtype(series.dtype)
            or pd.api.types.is_string_dtype(series.dtype)
        ):
            continue
        flags = [isinstance(v, str) and bool(_FORMULA_START.match(v)) for v in series.tolist()]
        if any(flags):
            mask = pd.Series(flags, index=series.index)
            out.isetitem(pos, series.where(~mask, "'" + series.astype("string")))
    return out


def write_csv(frame: pd.DataFrame, path: str | Path) -> None:
    """Write a DataFrame as UTF-8 CSV, creating parent directories.

    Text cells that would be executed as spreadsheet formulas are neutralised (see
    :func:`neutralise_formulas`).
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    neutralise_formulas(frame).to_csv(target, index=False, encoding="utf-8")


def write_json(payload: Any, path: str | Path) -> None:
    """Write deterministic, readable JSON."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8"
    )


def sha256_file(path: str | Path) -> str:
    """Return a streaming SHA-256 digest for a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def maybe_path_hash(path: str | Path | None) -> dict[str, str] | None:
    """Return input path and digest, or None for an omitted input."""
    if path is None:
        return None
    candidate = Path(path)
    if not candidate.exists():
        return {"path": str(candidate), "sha256": "MISSING"}
    return {"path": str(candidate), "sha256": sha256_file(candidate)}
