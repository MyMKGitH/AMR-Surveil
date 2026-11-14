"""EpiResNet adapters.

Predictions are preserved as model outputs and are never converted into observed AST.

What EpiResNet v5 actually writes (inspected at commit 4cc89fa)
---------------------------------------------------------------
``test_predictions.csv`` holds the **held-out test split** only: every input column plus derived
columns (``mic_mg_l``, ``phenotype``, ``binary_valid``, ...) and the model outputs
``probability_resistant`` (temperature-calibrated) and ``prediction_resistant`` (0/1 at the
validation-selected threshold). Both are ``NaN`` where no breakpoint-derived label exists. Rows
are **protein-antibiotic pairs**, not isolate-level phenotype calls. The training script does not
export predicted MIC; the Streamlit app computes it at inference only. ``final_metrics.json`` and
``run_manifest.json`` carry the threshold, temperature, split strategy and training-data hash.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pandas as pd

from amr_surveil.adapters.common import map_columns
from amr_surveil.constants import PREDICTION_COLUMNS
from amr_surveil.exceptions import ConfigurationError, InputFileError
from amr_surveil.io import read_csv, sha256_file, write_csv

PREDICTION_TARGETS = set(PREDICTION_COLUMNS)


def normalize_epiresnet_predictions(
    input_path: str | Path, output_path: str | Path, mapping_path: str | Path
) -> pd.DataFrame:
    """Normalize any EpiResNet-like inference CSV using a JSON column mapping."""
    return map_columns(
        input_path,
        output_path,
        mapping_path,
        required_targets={"sample_id", "antimicrobial"},
        allowed_targets=PREDICTION_TARGETS,
        output_columns=PREDICTION_COLUMNS,
        default_source_system="EpiResNet",
    )


def _read_json(path: str | Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigurationError(f"Could not read {path}: {exc}") from exc
    return payload if isinstance(payload, dict) else {}


def normalize_epiresnet_test_predictions(
    predictions_csv: str | Path,
    output_path: str | Path,
    *,
    final_metrics_path: str | Path | None = None,
    run_manifest_path: str | Path | None = None,
    checkpoint_path: str | Path | None = None,
    model_version: str = "v5",
) -> dict[str, object]:
    """Convert EpiResNet ``test_predictions.csv`` into the separate predictions stream.

    ``sample_id`` is taken from ``isolate_id``: it is the analyst's responsibility that this is the
    same identifier namespace as the WGS ``sample_id`` (unmatched IDs are reported by the audit).
    The checkpoint is identified by SHA-256 without loading it (no PyTorch needed).
    """
    path = Path(predictions_csv)
    frame = read_csv(path, table="epiresnet_test_predictions")
    need = {"isolate_id", "probability_resistant", "prediction_resistant"}
    missing = sorted(need - set(frame.columns))
    if missing:
        raise InputFileError(f"{path}: missing EpiResNet column(s): {', '.join(missing)}")
    drug_col = next((c for c in ("antibiotic_name", "antibiotic_id") if c in frame.columns), None)
    if drug_col is None:
        raise InputFileError(f"{path}: needs antibiotic_name or antibiotic_id")

    metrics = _read_json(final_metrics_path)
    manifest = _read_json(run_manifest_path)
    selection = metrics.get("selection", {})
    config = manifest.get("config", {})
    prob = pd.to_numeric(frame["probability_resistant"], errors="coerce")
    pred = pd.to_numeric(frame["prediction_resistant"], errors="coerce")
    valid = (
        frame["binary_valid"].astype("string").str.strip().str.lower().isin(["true", "1", "yes"])
        if "binary_valid" in frame.columns
        else prob.notna()
    )
    keep = valid.fillna(False).to_numpy() & prob.notna().to_numpy() & pred.notna().to_numpy()
    n_without = int((~keep).sum())

    out = pd.DataFrame(
        {
            "sample_id": frame["isolate_id"],
            "antimicrobial": frame[drug_col],
            "predicted_mic": pd.NA,
            "predicted_mic_unit": pd.NA,
            "predicted_class": pred.map(
                lambda v: "predicted_resistant" if v == 1 else "predicted_susceptible"
            ),
            "prediction_score": prob.round(6),
            "model_name": "EpiResNet",
            "model_version": model_version,
            "checkpoint_id": sha256_file(checkpoint_path) if checkpoint_path else pd.NA,
            "training_data_version": manifest.get("input_sha256", pd.NA),
            "split_strategy": config.get("split_strategy", pd.NA),
            "source_system": "EpiResNet",
            "prediction_unit": "protein-antibiotic pair",
            "protein_cluster": frame["protein_cluster"]
            if "protein_cluster" in frame.columns
            else pd.NA,
            "decision_threshold": selection.get("selected_threshold_on_validation", pd.NA),
            "probability_temperature": selection.get("probability_temperature", pd.NA),
            "evaluation_split": "test",
        }
    )
    if "predicted_log2_mic" in frame.columns:
        log2 = pd.to_numeric(frame["predicted_log2_mic"], errors="coerce")
        out["predicted_mic"] = log2.map(
            lambda v: round(2.0**v, 6) if pd.notna(v) and math.isfinite(v) else pd.NA
        )
        out["predicted_mic_unit"] = out["predicted_mic"].map(
            lambda v: "mg/L" if pd.notna(v) else pd.NA
        )
    out = out[keep].reindex(columns=PREDICTION_COLUMNS).reset_index(drop=True)
    write_csv(out, output_path)
    return {
        "n_input_rows": len(frame),
        "n_predictions_written": len(out),
        "n_rows_without_model_output": n_without,
        "checkpoint_hashed": checkpoint_path is not None,
        "threshold_recorded": "selected_threshold_on_validation" in selection,
        "split_strategy": config.get("split_strategy"),
        "note": "Rows are protein-antibiotic pairs from the held-out test split; no isolate-level phenotype is implied.",
    }
