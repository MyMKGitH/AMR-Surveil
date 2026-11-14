import json

import pandas as pd

from amr_surveil import __version__
from amr_surveil.cli import build_parser, main
from amr_surveil.synthetic import generate_synthetic_dataset
from tests.test_adapters_real_formats import tsv, v3_row


def test_parser_exposes_all_commands():
    sub = next(a for a in build_parser()._actions if a.dest == "command")
    assert {
        "run",
        "audit",
        "qa",
        "demo",
        "normalize-fishamr",
        "normalize-amrfinder",
        "normalize-epiresnet",
    } <= set(sub.choices)


def test_version_flag(capsys):
    try:
        main(["--version"])
    except SystemExit as exc:
        assert exc.code == 0
    assert __version__ in capsys.readouterr().out


def test_demo_command_produces_report(tmp_path, capsys):
    assert main(["demo", "--out", str(tmp_path), "--samples", "40"]) == 0
    assert (tmp_path / "results" / "report.md").is_file() and (
        tmp_path / "inputs" / "config.yaml"
    ).is_file()
    assert "genotype_phenotype_decision" in capsys.readouterr().out


def test_run_fails_with_status_2_on_error_findings_and_qa_status_1(tmp_path, capsys):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=30, inject_defects=True)
    args = ["--samples", str(paths["samples"]), "--ast", str(paths["ast"])]
    assert main(["run", *args, "--out", str(tmp_path / "r")]) == 2
    assert "Strict quality gate" in capsys.readouterr().err
    assert (
        main(["qa", *args, "--out", str(tmp_path / "q")]) == 1
    )  # CI-friendly: findings but outputs written
    assert (tmp_path / "q" / "audit_issues.csv").is_file()
    assert main(["audit", *args, "--out", str(tmp_path / "a")]) == 0


def test_normalize_amrfinder_command(tmp_path, capsys):
    rpt = tsv(tmp_path / "S1.tsv", [v3_row("S1", "blaTEM-1")])
    code = main(
        [
            "normalize-amrfinder",
            "--input",
            str(rpt),
            "--output",
            str(tmp_path / "c.csv"),
            "--screening-output",
            str(tmp_path / "s.csv"),
            "--tool-version",
            "3.12.8",
        ]
    )
    assert code == 0 and json.loads(capsys.readouterr().out)["n_calls_written"] == 1


def test_bad_input_gives_clean_error_not_traceback(tmp_path, capsys):
    code = main(
        [
            "normalize-amrfinder",
            "--input",
            str(tmp_path / "missing.tsv"),
            "--output",
            str(tmp_path / "c.csv"),
        ]
    )
    assert code == 2 and "ERROR" in capsys.readouterr().err
    assert main(["normalize-fishamr", "--input", "x", "--output", "y", "--kind", "amr"]) == 2
    assert main(["run", "--samples", str(tmp_path / "nope.csv"), "--out", str(tmp_path / "o")]) == 2


def test_epiresnet_cli_requires_mapping_or_flag(tmp_path, capsys):
    f = tmp_path / "x.csv"
    pd.DataFrame({"a": [1]}).to_csv(f, index=False)
    assert (
        main(["normalize-epiresnet", "--input", str(f), "--output", str(tmp_path / "o.csv")]) == 2
    )
