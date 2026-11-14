"""UI-independent helpers used by ``app.py`` (kept out of the Streamlit script so they can be tested)."""

from __future__ import annotations

import io
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from amr_surveil.cli import DEMO_CONFIG
from amr_surveil.config import validate_config
from amr_surveil.constants import DEDUPLICATION_RULES, EVIDENCE_TIERS
from amr_surveil.datasources.common import read_table_upload as read_table_upload
from amr_surveil.pipeline import run_pipeline
from amr_surveil.synthetic import generate_synthetic_dataset

#: (key, upload label, required)
UPLOAD_SLOTS: list[tuple[str, str, bool]] = [
    ("samples", "samples.csv", True),
    ("amr_calls", "amr_calls.csv", False),
    ("screening", "screening.csv (which samples were genotyped)", False),
    ("ast", "ast.csv", False),
    ("breakpoints", "breakpoints.csv", False),
    ("gene_drug_map", "gene_drug_map.csv (cited)", False),
    ("aggregate", "aggregate_surveillance.csv", False),
    ("predictions", "model_predictions.csv", False),
    ("plasmids", "plasmid_calls.csv", False),
    ("incidence", "incidence.csv", False),
    ("consumption", "consumption.csv", False),
]
_ARG = {
    "samples": "samples_path",
    "amr_calls": "amr_calls_path",
    "screening": "screening_path",
    "ast": "ast_path",
    "breakpoints": "breakpoints_path",
    "aggregate": "aggregate_path",
    "predictions": "predictions_path",
    "plasmids": "plasmids_path",
    "incidence": "incidence_path",
    "consumption": "consumption_path",
    "gene_drug_map": "gene_drug_map_path",
}


def build_config(
    *,
    strict: bool = True,
    deduplication_rule: str = "none",
    excluded_evidence_tiers: list[str] | None = None,
    study_name: str = "",
    sample_selection_documented: bool = False,
    selection_description: str = "",
    predictions_independence_attested: bool = False,
) -> dict[str, dict[str, Any]]:
    """Validated configuration from UI choices (unknown values raise ConfigurationError)."""
    tiers = (
        list(EVIDENCE_TIERS[1:3]) if excluded_evidence_tiers is None else excluded_evidence_tiers
    )
    return validate_config(
        {
            "analysis": {
                "deduplication_rule": deduplication_rule,
                "excluded_evidence_tiers": tiers,
            },
            "quality": {"strict": strict},
            "study": {
                "name": study_name,
                "sample_selection_documented": sample_selection_documented,
                "selection_description": selection_description,
                "predictions_independence_attested": predictions_independence_attested,
            },
        }
    )


def run_from_files(
    paths: dict[str, Path], config: dict[str, dict[str, Any]], out_dir: Path
) -> dict[str, Any]:
    """Run the pipeline on uploaded files with a UI-built configuration."""
    out_dir.mkdir(parents=True, exist_ok=True)
    config_path = out_dir.parent / "ui_config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    kwargs = {_ARG[k]: v for k, v in paths.items() if k in _ARG}
    return run_pipeline(**kwargs, config_path=config_path, output_dir=out_dir)


def run_demo(workdir: Path, n_samples: int = 120) -> tuple[dict[str, Any], Path]:
    inputs = generate_synthetic_dataset(workdir / "inputs", n_samples=n_samples)
    config_path = workdir / "demo_config.yaml"
    config_path.write_text(yaml.safe_dump(DEMO_CONFIG, sort_keys=False), encoding="utf-8")
    kwargs = {_ARG[k]: v for k, v in inputs.items() if k in _ARG}
    out = workdir / "results"
    return run_pipeline(**kwargs, config_path=config_path, output_dir=out), out


def read_result(result_dir: Path, name: str) -> pd.DataFrame:
    """Read an output CSV; a missing or header-less file is an empty table, never a crash."""
    path = result_dir / name
    if not path.is_file():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def zip_outputs(folder: Path) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(folder.rglob("*")):
            if path.is_file():
                archive.write(path, arcname=path.relative_to(folder))
    return buffer.getvalue()


def error_bars(frame: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Upper/lower Plotly error-bar lengths from Wilson bounds (never negative)."""
    up = (frame["wilson_high_pct"] - frame["resistant_pct"]).clip(lower=0)
    down = (frame["resistant_pct"] - frame["wilson_low_pct"]).clip(lower=0)
    return up, down


DEDUP_CHOICES = list(DEDUPLICATION_RULES)
TIER_CHOICES = list(EVIDENCE_TIERS)


# --------------------------------------------------------------------------- public data hub
def analyse_public_dataset(
    dataset: Any, workdir: Path, config: dict[str, dict[str, Any]]
) -> tuple[dict[str, Any], Path]:
    """Write a public dataset as canonical inputs and run the standard pipeline on it."""
    from amr_surveil.datasources import write_inputs

    paths = write_inputs(dataset, workdir / "inputs")
    if "samples" not in paths or "ast" not in paths:
        raise ValueError("This dataset has no isolate-level samples and AST records to analyse.")
    out = workdir / "results"
    return run_from_files({k: v for k, v in paths.items() if k in _ARG}, config, out), out


def public_inputs_zip(dataset: Any, workdir: Path) -> bytes:
    from amr_surveil.datasources import write_inputs

    write_inputs(dataset, workdir / "inputs")
    return zip_outputs(workdir / "inputs")


def context_labels(frame: pd.DataFrame) -> list[str]:
    """Distinct 'method | standard version' labels present in an indicator table."""
    if frame.empty:
        return []
    cols = [
        c for c in ("method", "breakpoint_standard", "breakpoint_version") if c in frame.columns
    ]
    if not cols:
        return []
    labels = frame[cols].fillna("").astype(str).agg(" | ".join, axis=1)
    return sorted(labels.unique())


def filter_context(frame: pd.DataFrame, label: str | None) -> pd.DataFrame:
    if frame.empty or not label:
        return frame
    cols = [
        c for c in ("method", "breakpoint_standard", "breakpoint_version") if c in frame.columns
    ]
    labels = frame[cols].fillna("").astype(str).agg(" | ".join, axis=1)
    return frame[labels == label]
