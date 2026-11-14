from pathlib import Path

import pandas as pd
import pytest

from amr_surveil.exceptions import QualityGateError
from amr_surveil.pipeline import run_pipeline
from amr_surveil.synthetic import generate_synthetic_dataset


def test_synthetic_end_to_end_outputs_are_complete(tmp_path: Path):
    paths = generate_synthetic_dataset(tmp_path / "inputs", n_samples=12, seed=7)
    out = tmp_path / "results"
    summary = run_pipeline(
        samples_path=paths["samples"],
        amr_calls_path=paths["amr_calls"],
        ast_path=paths["ast"],
        aggregate_path=paths["aggregate"],
        predictions_path=paths["predictions"],
        plasmids_path=paths["plasmids"],
        incidence_path=paths["incidence"],
        consumption_path=paths["consumption"],
        output_dir=out,
    )
    expected = {
        "audit_issues.csv",
        "resistance_by_group.csv",
        "aggregate_surveillance_context.csv",
        "incidence_context.csv",
        "consumption_context.csv",
        "matched_isolates.csv",
        "model_prediction_context.csv",
        "plasmid_calls_linked.csv",
        "gene_drug_descriptive_stratification.csv",
        "analysis_summary.json",
        "run_manifest.json",
        "report.md",
        "exclusion_log.csv",
        "data_inventory.csv",
        "readiness_assessment.json",
        "comparability_flags.csv",
    }
    assert expected.issubset({p.name for p in out.iterdir()})
    assert summary["n_samples"] == 12
    assert summary["n_predictions"] > 0
    assert summary["n_aggregate_rows"] == 3
    assert summary["n_incidence_rows"] == 1
    assert summary["n_consumption_rows"] == 1
    assert "Separate surveillance streams" in (out / "report.md").read_text()
    assert (out / "report.md").read_text().startswith("# AMR-Surveil analysis report")


def test_aggregate_remains_separate(tmp_path: Path):
    paths = generate_synthetic_dataset(tmp_path / "inputs", n_samples=8)
    out = tmp_path / "results"
    run_pipeline(
        samples_path=paths["samples"],
        amr_calls_path=paths["amr_calls"],
        ast_path=paths["ast"],
        aggregate_path=paths["aggregate"],
        output_dir=out,
    )
    aggregate = pd.read_csv(out / "aggregate_surveillance_context.csv")
    assert set(aggregate["evidence_level"]) == {"aggregate surveillance; not isolate-linked"}
    assert "sample_id" not in aggregate.columns


def test_strict_gate_stops_invalid_ast_and_writes_audit(tmp_path: Path):
    paths = generate_synthetic_dataset(tmp_path / "inputs", n_samples=8)
    ast = pd.read_csv(paths["ast"])
    ast.loc[0, "result"] = "Intermediate"
    bad_ast = tmp_path / "bad_ast.csv"
    ast.to_csv(bad_ast, index=False)
    out = tmp_path / "results"
    with pytest.raises(QualityGateError):
        run_pipeline(samples_path=paths["samples"], ast_path=bad_ast, output_dir=out)
    issues = pd.read_csv(out / "audit_issues.csv")
    assert "invalid_ast_category" in set(issues["code"])


def test_duplicate_sample_ids_fail_strict_gate(tmp_path: Path):
    samples_path = tmp_path / "samples.csv"
    pd.DataFrame({"sample_id": ["A", "A"], "organism": ["E. coli", "E. coli"]}).to_csv(
        samples_path, index=False
    )
    with pytest.raises(QualityGateError):
        run_pipeline(samples_path=samples_path, output_dir=tmp_path / "results")


def test_fractional_aggregate_counts_fail_strict_gate(tmp_path: Path):
    paths = generate_synthetic_dataset(tmp_path / "inputs", n_samples=8)
    aggregate = pd.read_csv(paths["aggregate"])
    aggregate["resistant_n"] = aggregate["resistant_n"].astype(float)
    aggregate.loc[0, "resistant_n"] = 1.5
    bad = tmp_path / "aggregate_fractional.csv"
    aggregate.to_csv(bad, index=False)
    with pytest.raises(QualityGateError):
        run_pipeline(
            samples_path=paths["samples"],
            aggregate_path=bad,
            output_dir=tmp_path / "results",
        )
    issues = pd.read_csv(tmp_path / "results" / "audit_issues.csv")
    assert "non_integer_aggregate_counts" in set(issues["code"])


def test_invalid_boolean_config_is_rejected_instead_of_truthiness_coercion(tmp_path: Path):
    from amr_surveil.config import load_config as _load_config
    from amr_surveil.exceptions import ConfigurationError

    config = tmp_path / "bad.yaml"
    config.write_text("quality:\n  strict: 'false'\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="quality.strict must be a YAML boolean"):
        _load_config(config)


def test_unknown_config_keys_are_rejected(tmp_path: Path):
    from amr_surveil.config import load_config as _load_config
    from amr_surveil.exceptions import ConfigurationError

    config = tmp_path / "bad.yaml"
    config.write_text("analysis:\n  confidence: 0.95\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="Unsupported analysis setting"):
        _load_config(config)


def test_invalid_population_denominator_stops_pipeline(tmp_path: Path):
    paths = generate_synthetic_dataset(tmp_path / "inputs", n_samples=8)
    incidence = pd.read_csv(paths["incidence"])
    incidence.loc[0, "population_n"] = 0
    bad = tmp_path / "bad_incidence.csv"
    incidence.to_csv(bad, index=False)
    with pytest.raises(QualityGateError):
        run_pipeline(
            samples_path=paths["samples"],
            incidence_path=bad,
            output_dir=tmp_path / "results",
        )
    issues = pd.read_csv(tmp_path / "results" / "audit_issues.csv")
    assert "zero_population_denominator" in set(issues["code"])


def test_invalid_consumption_period_stops_pipeline(tmp_path: Path):
    paths = generate_synthetic_dataset(tmp_path / "inputs", n_samples=8)
    consumption = pd.read_csv(paths["consumption"])
    consumption.loc[0, "period_days"] = 0
    bad = tmp_path / "bad_consumption.csv"
    consumption.to_csv(bad, index=False)
    with pytest.raises(QualityGateError):
        run_pipeline(
            samples_path=paths["samples"],
            consumption_path=bad,
            output_dir=tmp_path / "results",
        )
    issues = pd.read_csv(tmp_path / "results" / "audit_issues.csv")
    assert "zero_period_days" in set(issues["code"])


def test_audit_mode_preserves_effective_config_in_manifest(tmp_path: Path):
    import json

    from amr_surveil.cli import main

    paths = generate_synthetic_dataset(tmp_path / "inputs", n_samples=8)
    bad_ast = pd.read_csv(paths["ast"])
    bad_ast.loc[0, "result"] = "Intermediate"
    bad_path = tmp_path / "bad_ast.csv"
    bad_ast.to_csv(bad_path, index=False)
    out = tmp_path / "exploratory"
    status = main(
        ["audit", "--samples", str(paths["samples"]), "--ast", str(bad_path), "--out", str(out)]
    )
    assert status == 0
    assert not (out / "exploratory_config.yaml").exists()  # no temporary config files any more
    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["effective_configuration"]["quality"]["strict"] is False
    assert "audit_issues.csv" in manifest["outputs_sha256"]
