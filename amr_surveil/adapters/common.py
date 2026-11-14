"""Shared adapter mapping utilities."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from amr_surveil.exceptions import ConfigurationError
from amr_surveil.io import read_csv, write_csv


def load_mapping(path: str | Path) -> dict[str, Any]:
    """Load and validate a simple JSON column-mapping file."""
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigurationError(f"Could not load adapter mapping {path}: {exc}") from exc
    if not isinstance(payload, dict) or not isinstance(payload.get("columns"), dict):
        raise ConfigurationError("Adapter mapping must contain a 'columns' object")
    payload.setdefault("defaults", {})
    if not isinstance(payload["defaults"], dict):
        raise ConfigurationError("Adapter mapping 'defaults' must be an object")
    return payload


def map_columns(
    input_path: str | Path,
    output_path: str | Path,
    mapping_path: str | Path,
    *,
    required_targets: set[str],
    allowed_targets: set[str],
    output_columns: list[str],
    default_source_system: str,
) -> pd.DataFrame:
    """Map arbitrary exported column names to a canonical output contract."""
    payload = load_mapping(mapping_path)
    columns = payload["columns"]
    unsupported = set(columns) - allowed_targets
    missing = required_targets - set(columns)
    if unsupported:
        raise ConfigurationError(f"Unsupported target field(s): {', '.join(sorted(unsupported))}")
    if missing:
        raise ConfigurationError(
            f"Mapping missing required target field(s): {', '.join(sorted(missing))}"
        )
    source = read_csv(input_path, table="adapter_input")
    output = pd.DataFrame(index=source.index)
    for target, source_col in columns.items():
        if source_col not in source.columns:
            raise ConfigurationError(
                f"Mapped source column '{source_col}' for '{target}' is absent from input"
            )
        output[target] = source[source_col]
    defaults = payload["defaults"]
    for col, value in defaults.items():
        if col not in allowed_targets:
            raise ConfigurationError(f"Unsupported default field: {col}")
        if col not in output.columns:
            output[col] = value
    if "source_system" in allowed_targets and "source_system" not in output.columns:
        output["source_system"] = default_source_system
    if "model_name" in allowed_targets and "model_name" not in output.columns:
        output["model_name"] = default_source_system
    for col in output_columns:
        if col not in output.columns:
            output[col] = pd.NA
    output = output[output_columns]
    write_csv(output, output_path)
    return output
