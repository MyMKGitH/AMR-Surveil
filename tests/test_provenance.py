"""verify and templates commands."""

import json

import pandas as pd
import pytest

from amr_surveil import constants as C
from amr_surveil.cli import main
from amr_surveil.config import load_config
from amr_surveil.exceptions import ConfigurationError, InputFileError
from amr_surveil.pipeline import run_pipeline
from amr_surveil.provenance import TEMPLATE_TABLES, verify_run, write_templates
from tests.conftest import pipeline_kwargs


def _run(demo, out):
    return run_pipeline(
        **pipeline_kwargs(demo["paths"], config_path=demo["config"], output_dir=out)
    )


def test_verify_passes_on_fresh_run_and_detects_every_kind_of_tampering(tmp_path, demo):
    out = tmp_path / "res"
    _run(demo, out)
    assert verify_run(out)["ok"] is True

    target = out / "resistance_by_group.csv"
    original = target.read_text(encoding="utf-8")
    target.write_text(original.replace("0.0", "0.1", 1), encoding="utf-8")
    report = verify_run(out)
    assert report["ok"] is False and report["mismatched"] == ["resistance_by_group.csv"]
    target.write_text(original, encoding="utf-8")
    assert verify_run(out)["ok"] is True

    (out / "report.md").unlink()
    assert verify_run(out)["missing"] == ["report.md"]
    (out / "report.md").write_text("x", encoding="utf-8")
    (out / "extra.txt").write_text("x", encoding="utf-8")
    assert "extra.txt" in verify_run(out)["unlisted"]


def test_verify_checks_inputs_when_asked(tmp_path, demo):
    out = tmp_path / "res"
    _run(demo, out)
    report = verify_run(out, check_inputs=True)
    assert report["ok"] and report["inputs"]["samples"] == "ok"


def test_verify_rejects_missing_or_corrupt_manifest(tmp_path):
    with pytest.raises(InputFileError, match="No run_manifest"):
        verify_run(tmp_path)
    (tmp_path / "run_manifest.json").write_text("{not json", encoding="utf-8")
    with pytest.raises(InputFileError, match="not a valid"):
        verify_run(tmp_path)
    (tmp_path / "run_manifest.json").write_text(json.dumps({"x": 1}), encoding="utf-8")
    with pytest.raises(InputFileError, match="not a valid"):
        verify_run(tmp_path)


def test_manifest_records_schema_version(tmp_path, demo):
    out = tmp_path / "res"
    _run(demo, out)
    assert json.loads((out / "run_manifest.json").read_text())["manifest_schema_version"] == 1


def test_cli_verify_exit_codes(tmp_path, capsys):
    out = tmp_path / "demo"
    assert main(["demo", "--out", str(out), "--samples", "30"]) == 0
    assert main(["verify", str(out / "results")]) == 0
    (out / "results" / "report.md").write_text("tampered", encoding="utf-8")
    assert main(["verify", str(out / "results")]) == 1
    assert main(["verify", str(tmp_path / "nowhere")]) == 2


def test_templates_match_contracts_and_config_validates(tmp_path):
    written = write_templates(tmp_path / "t")
    assert len(written) == len(TEMPLATE_TABLES) + 1
    for name, columns in TEMPLATE_TABLES.items():
        assert pd.read_csv(tmp_path / "t" / name).columns.tolist() == columns
    cfg = load_config(tmp_path / "t" / "config.yaml")
    assert cfg["quality"]["strict"] is True and cfg["analysis"]["deduplication_rule"] == "none"
    assert set(TEMPLATE_TABLES["samples.csv"]) == set(C.SAMPLE_COLUMNS)


def test_templates_never_overwrite_silently(tmp_path):
    write_templates(tmp_path)
    (tmp_path / "samples.csv").write_text("sample_id\nKEEP\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="Refusing to overwrite"):
        write_templates(tmp_path)
    assert "KEEP" in (tmp_path / "samples.csv").read_text()
    write_templates(tmp_path, force=True)
    assert "KEEP" not in (tmp_path / "samples.csv").read_text()


def test_empty_templates_do_not_crash_the_pipeline(tmp_path):
    t = tmp_path / "t"
    write_templates(t)
    out = tmp_path / "o"
    code = main(
        [
            "audit",
            "--samples",
            str(t / "samples.csv"),
            "--ast",
            str(t / "ast.csv"),
            "--out",
            str(out),
        ]
    )
    assert code == 0 and (out / "run_manifest.json").exists()
