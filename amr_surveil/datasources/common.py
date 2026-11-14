"""Shared helpers: result vocabularies, sector inference, and the :class:`PublicDataset` container.

Everything converting public records into AMR-Surveil inputs follows the same conservative rules:
a value that cannot be mapped is **dropped and counted** (never guessed), and the loss is reported in
``PublicDataset.notes`` so it reaches the user and the saved provenance.
"""

from __future__ import annotations

import io
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from amr_surveil.constants import AST_COLUMNS, SAMPLE_COLUMNS

#: Public-source phenotype wording -> S/I/R. ``nonsusceptible`` is deliberately absent: it merges I and
#: R and so cannot be assigned to either without inventing information.
RESULT_MAP = {
    "r": "R",
    "resistant": "R",
    "s": "S",
    "susceptible": "S",
    "sensitive": "S",
    "i": "I",
    "intermediate": "I",
    "susceptible increased exposure": "I",
    "susceptible dose dependent": "I",
    "sdd": "I",
    "intermediate resistance": "I",
}

_HUMAN = re.compile(r"(?<![a-z-])(homo[ _]sapiens|h\. ?sapiens|humans?|patients?)\b", re.I)
_NONHUMAN = re.compile(r"\bnon[- ]?human\b", re.I)
_ENVIRONMENT = re.compile(
    r"\b(environment|environmental|soil|sediment|sewage|wastewater|water|surface|biofilm)\b", re.I
)


def _clean(value: object) -> str:
    """Blank for missing/NaN/placeholder values, otherwise the stripped string."""
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {
        "",
        "nan",
        "none",
        "null",
        "na",
        "n/a",
        "missing",
        "not available",
        "not collected",
        "not applicable",
        "unknown",
        "not provided",
        "not determined",
        "not recorded",
        "restricted access",
        "unspecified",
        "-",
        "--",
        "?",
    }:
        return ""
    return text


def clean(value: object) -> str:
    return _clean(value)


def map_result(value: object) -> str | None:
    text = re.sub(r"[^a-z ]+", " ", _clean(value).lower())
    return RESULT_MAP.get(re.sub(r"\s+", " ", text).strip())


def infer_sector(host: object, source: object = "", epi_type: object = "") -> str:
    """Coarse sector from host / isolation source. Ambiguity returns ``unknown`` (never a guess)."""
    host_text, source_text = _clean(host), _clean(source)
    if host_text and not _NONHUMAN.search(host_text) and _HUMAN.search(host_text):
        return "human"
    if (
        host_text
        and _ENVIRONMENT.search(host_text)
        and not re.search(r"\bsapiens\b", host_text, re.I)
    ):
        return "environment"
    if host_text:
        return "animal"
    if _HUMAN.search(source_text):
        return "human"
    if _ENVIRONMENT.search(source_text) or _clean(epi_type).lower().startswith("environmental"):
        return "environment"
    return "unknown"


def normalise_date(value: object) -> str:
    """ISO date prefix (YYYY, YYYY-MM or YYYY-MM-DD) from a free-form collection date, else blank."""
    text = _clean(value)
    match = re.match(r"^(\d{4})(?:-(\d{2})(?:-(\d{2}))?)?", text)
    if not match:
        return ""
    year, month, day = match.groups()
    if month and not 1 <= int(month) <= 12:
        return year
    if day and not 1 <= int(day) <= 31:
        return f"{year}-{month}"
    return "-".join(p for p in (year, month, day) if p)


def title_drug(name: object) -> str:
    text = _clean(name)
    return text[:1].upper() + text[1:].lower() if text else ""


@dataclass
class PublicDataset:
    """AMR-Surveil inputs derived from one public source, with the provenance to cite it."""

    source_id: str
    samples: pd.DataFrame
    ast: pd.DataFrame
    provenance: dict[str, Any]
    notes: list[str] = field(default_factory=list)
    aggregate: pd.DataFrame | None = None
    amr_calls: pd.DataFrame | None = None
    screening: pd.DataFrame | None = None

    def summary(self) -> dict[str, Any]:
        return {
            "source": self.source_id,
            "n_samples": int(len(self.samples)),
            "n_ast_records": int(len(self.ast)),
            "n_aggregate_rows": 0 if self.aggregate is None else int(len(self.aggregate)),
            "notes": list(self.notes),
        }


def write_inputs(dataset: PublicDataset, out_dir: str | Path) -> dict[str, Path]:
    """Write the dataset as canonical AMR-Surveil CSVs plus ``public_data_provenance.json``.

    The returned mapping uses the keys of ``ui_support.UPLOAD_SLOTS`` so it can feed
    ``run_from_files`` / the CLI unchanged.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    tables = {
        "samples": dataset.samples,
        "ast": dataset.ast,
        "aggregate": dataset.aggregate,
        "amr_calls": dataset.amr_calls,
        "screening": dataset.screening,
    }
    for key, frame in tables.items():
        if frame is None or frame.empty:
            continue
        path = out / f"{key}.csv"
        frame.to_csv(path, index=False)
        paths[key] = path
    record = {**dataset.provenance, "notes": dataset.notes, "summary": dataset.summary()}
    (out / "public_data_provenance.json").write_text(
        json.dumps(record, indent=2, default=str), encoding="utf-8"
    )
    return paths


def empty_samples() -> pd.DataFrame:
    return pd.DataFrame(columns=SAMPLE_COLUMNS)


def empty_ast() -> pd.DataFrame:
    return pd.DataFrame(columns=AST_COLUMNS)


def header_key(name: object) -> str:
    """Lower-case alphanumeric header key so 'Genome ID', 'genome_id' and 'GenomeId' all match."""
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def pick(frame: pd.DataFrame, *aliases: str) -> pd.Series:
    """First column whose normalised header matches an alias; an all-blank series if none do."""
    lookup = {header_key(c): c for c in frame.columns}
    for alias in aliases:
        if header_key(alias) in lookup:
            return frame[lookup[header_key(alias)]]
    return pd.Series([""] * len(frame), index=frame.index, dtype="object")


def read_table_upload(data: bytes, filename: str = "") -> pd.DataFrame:
    """Parse an uploaded CSV/TSV (delimiter sniffed from the header) with every cell kept as text."""
    text = data.decode("utf-8-sig", errors="replace")
    if not text.strip():
        return pd.DataFrame()
    header = text.splitlines()[0] if text.strip() else ""
    sep = (
        "\t"
        if filename.lower().endswith((".tsv", ".txt")) or header.count("\t") > header.count(",")
        else ","
    )
    return pd.read_csv(io.StringIO(text), sep=sep, dtype=str, keep_default_na=False)
