"""Regression tests for defects found by adversarial probing of v0.1.0 (see docs/VALIDATION.md)."""

import pandas as pd
import pytest

from amr_surveil.pipeline import run_pipeline

MANIFEST = (
    "sample_id,organism,collection_date,country,source_sector,patient_id\n"
    "A1,Escherichia coli,2024-01-05,NO,human,P1\n"
    "A2,Escherichia coli,05/06/2024,NO,human,P2\n"
    "A4,Escherichia coli,2024-02-01,NO,human,P4\n"
)
AST_HEAD = "sample_id,antimicrobial,result,mic,method,breakpoint_standard,breakpoint_version\n"


def run(tmp_path, ast=None, aggregate=None, manifest=MANIFEST):
    (tmp_path / "s.csv").write_text(manifest)
    kw = {}
    if ast:
        (tmp_path / "a.csv").write_text(AST_HEAD + ast)
        kw["ast_path"] = tmp_path / "a.csv"
    if aggregate:
        (tmp_path / "g.csv").write_text(aggregate)
        kw["aggregate_path"] = tmp_path / "g.csv"
    summary = run_pipeline(samples_path=tmp_path / "s.csv", output_dir=tmp_path / "out", **kw)
    return (
        summary,
        pd.read_csv(tmp_path / "out" / "audit_issues.csv"),
        pd.read_csv(tmp_path / "out" / "resistance_by_group.csv"),
    )


def test_A_blank_result_is_a_warning_not_a_gate_stopping_error(tmp_path):
    _, issues, ind = run(
        tmp_path, ast="A1,Cefotaxime,,1,BMD,EUCAST,16.1\nA4,Cefotaxime,R,4,BMD,EUCAST,16.1\n"
    )
    row = issues[issues.code == "missing_result"].iloc[0]
    assert row.severity == "WARNING" and "invalid_ast_category" not in set(issues.code)
    assert ind.n_tested.sum() == 1  # the blank record is excluded from the denominator, and logged


def test_B_censored_mic_is_understood_not_flagged_invalid(tmp_path):
    _, issues, _ = run(
        tmp_path, ast="A1,Cefotaxime,S,<=0.5,BMD,EUCAST,16.1\nA4,Cefotaxime,R,>32,BMD,EUCAST,16.1\n"
    )
    assert "invalid_mic" not in set(issues.code) and "censored_mic_present" in set(issues.code)


def test_C_label_case_variants_are_merged_and_reported(tmp_path):
    _, issues, ind = run(
        tmp_path,
        ast="A1,Cefotaxime,R,,BMD,EUCAST,16.1\nA2,cefotaxime,R,,BMD,EUCAST,16.1\nA4,CEFOTAXIME ,S,,BMD,EUCAST,16.1\n",
    )
    assert "inconsistent_antimicrobial_labels" in set(issues.code)
    assert ind.antimicrobial.nunique() == 1 and ind.n_tested.sum() == 3  # one stratum, not three


def test_D_duplicate_aggregate_rows_stop_the_strict_gate(tmp_path):
    from amr_surveil.exceptions import QualityGateError

    row = "2024,NO,Escherichia coli,Cefotaxime,10,100,SRC,def\n"
    head = (
        "year,country,organism,antimicrobial,resistant_n,tested_n,source_name,metric_definition\n"
    )
    (tmp_path / "s.csv").write_text(MANIFEST)
    (tmp_path / "g.csv").write_text(head + row + row)
    with pytest.raises(QualityGateError, match="duplicate_aggregate_key"):
        run_pipeline(
            samples_path=tmp_path / "s.csv",
            aggregate_path=tmp_path / "g.csv",
            output_dir=tmp_path / "o",
        )


def test_E_future_and_ambiguous_dates_are_reported(tmp_path):
    manifest = MANIFEST + "A5,Escherichia coli,2099-01-01,NO,human,P5\n"
    (tmp_path / "s.csv").write_text(manifest)
    s = run_pipeline(
        samples_path=tmp_path / "s.csv",
        output_dir=tmp_path / "o",
        config_overrides={"quality": {"strict": False}},
    )
    issues = pd.read_csv(tmp_path / "o" / "audit_issues.csv")
    assert {"future_collection_date", "non_iso_collection_date"} <= set(issues.code)
    assert s["n_qa_errors"] == 1
