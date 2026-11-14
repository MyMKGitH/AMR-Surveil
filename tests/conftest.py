from pathlib import Path

import pytest
import yaml

from amr_surveil.cli import DEMO_CONFIG
from amr_surveil.pipeline import run_pipeline
from amr_surveil.synthetic import generate_synthetic_dataset


def pipeline_kwargs(paths: dict[str, Path], **extra):
    mapping = {
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
    kwargs = {arg: paths[key] for key, arg in mapping.items() if key in paths}
    kwargs.update(extra)
    return kwargs


@pytest.fixture(scope="session")
def demo(tmp_path_factory):
    """A full synthetic run (n=120) shared by read-only tests."""
    root = tmp_path_factory.mktemp("demo")
    paths = generate_synthetic_dataset(root / "inputs", n_samples=120)
    config = root / "config.yaml"
    config.write_text(yaml.safe_dump(DEMO_CONFIG), encoding="utf-8")
    summary = run_pipeline(
        **pipeline_kwargs(paths, config_path=config, output_dir=root / "results")
    )
    return {"paths": paths, "config": config, "out": root / "results", "summary": summary}
