"""Re-verify a finished run against its ``run_manifest.json`` and scaffold input templates."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml

from amr_surveil import constants as C
from amr_surveil.config import DEFAULTS
from amr_surveil.exceptions import ConfigurationError, InputFileError
from amr_surveil.io import sha256_file

TEMPLATE_TABLES: dict[str, list[str]] = {
    "samples.csv": C.SAMPLE_COLUMNS,
    "amr_calls.csv": C.AMR_CALL_COLUMNS,
    "screening.csv": C.SCREENING_COLUMNS,
    "ast.csv": C.AST_COLUMNS,
    "breakpoints.csv": C.BREAKPOINT_COLUMNS,
    "aggregate_surveillance.csv": C.AGGREGATE_COLUMNS,
    "model_predictions.csv": C.PREDICTION_COLUMNS,
    "plasmid_calls.csv": C.PLASMID_COLUMNS,
    "incidence.csv": C.INCIDENCE_COLUMNS,
    "consumption.csv": C.CONSUMPTION_COLUMNS,
    "gene_drug_map.csv": C.GENE_MAP_COLUMNS,
}


def verify_run(results_dir: str | Path, *, check_inputs: bool = False) -> dict[str, Any]:
    """Compare every recorded output hash (and optionally input hash) with the files on disk.

    Returns a dictionary with ``ok`` and the lists ``mismatched``, ``missing`` and ``unlisted``.
    ``run_manifest.json`` itself cannot contain its own hash and is therefore not checked.
    """
    root = Path(results_dir)
    manifest_path = root / "run_manifest.json"
    if not manifest_path.is_file():
        raise InputFileError(f"No run_manifest.json in {root}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        recorded = dict(manifest["outputs_sha256"])
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise InputFileError(f"{manifest_path} is not a valid AMR-Surveil manifest: {exc}") from exc
    mismatched, missing = [], []
    for name, digest in sorted(recorded.items()):
        target = root / name
        if not target.is_file():
            missing.append(name)
        elif sha256_file(target) != digest:
            mismatched.append(name)
    unlisted = sorted(
        p.name
        for p in root.iterdir()
        if p.is_file() and p.name not in {*recorded, "run_manifest.json"}
    )
    inputs_report: dict[str, str] = {}
    if check_inputs:
        for key, item in sorted((manifest.get("inputs") or {}).items()):
            if not item:
                continue
            path = Path(item["path"])
            if not path.is_file():
                inputs_report[key] = "not_found"
            else:
                inputs_report[key] = "ok" if sha256_file(path) == item["sha256"] else "changed"
    ok = not (mismatched or missing) and all(v == "ok" for v in inputs_report.values())
    return {
        "ok": ok,
        "software_version": manifest.get("software_version"),
        "n_outputs_checked": len(recorded),
        "mismatched": mismatched,
        "missing": missing,
        "unlisted": unlisted,
        "inputs": inputs_report,
    }


def write_templates(out_dir: str | Path, *, force: bool = False) -> list[Path]:
    """Write header-only input templates and a default ``config.yaml`` (never overwrites silently)."""
    root = Path(out_dir)
    root.mkdir(parents=True, exist_ok=True)
    targets = [*TEMPLATE_TABLES, "config.yaml"]
    existing = [name for name in targets if (root / name).exists()]
    if existing and not force:
        raise ConfigurationError(
            f"Refusing to overwrite existing file(s) in {root}: {', '.join(existing)} (use --force)"
        )
    written = []
    for name, columns in TEMPLATE_TABLES.items():
        (root / name).write_text(",".join(columns) + "\n", encoding="utf-8")
        written.append(root / name)
    config = {section: dict(values) for section, values in DEFAULTS.items()}
    (root / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    written.append(root / "config.yaml")
    return written
