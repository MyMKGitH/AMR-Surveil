"""Regression tests for the second independent review round (see CHANGELOG 1.0.0)."""

import pandas as pd

from amr_surveil.adapters.amrfinder import normalize_amrfinder_reports
from amr_surveil.adapters.epiresnet import normalize_epiresnet_test_predictions
from amr_surveil.adapters.fishamr import normalize_fishamr_feature_table
from amr_surveil.pipeline import run_pipeline
from amr_surveil.synthetic import generate_synthetic_dataset
from tests.conftest import pipeline_kwargs
from tests.test_independent_oracle import CONFIG

SAMPLES = (
    "sample_id,organism,collection_date,country,source_sector\n"
    "A,Escherichia coli,2024-05-20,NO,human\n"
)


def run(tmp_path, samples=SAMPLES, ast=None, report=False, overrides=None, **tables):
    (tmp_path / "s.csv").write_text(samples, encoding="utf-8")
    kwargs = {}
    if ast is not None:
        (tmp_path / "a.csv").write_text(ast, encoding="utf-8")
        kwargs["ast_path"] = tmp_path / "a.csv"
    for key, text in tables.items():
        (tmp_path / f"{key}.csv").write_text(text, encoding="utf-8")
        kwargs[f"{key}_path"] = tmp_path / f"{key}.csv"
    summary = run_pipeline(
        samples_path=tmp_path / "s.csv",
        output_dir=tmp_path / "out",
        config_overrides={"quality": {"strict": False}, **(overrides or {})},
        **kwargs,
    )
    out = tmp_path / "out"
    issues = pd.read_csv(out / "audit_issues.csv")
    return out, summary, set(issues["code"])


def test_upper_case_evidence_tiers_are_still_excluded(tmp_path):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=60, seed=11)
    base = run_pipeline(
        **pipeline_kwargs(paths, output_dir=tmp_path / "b", config_overrides=CONFIG)
    )
    calls = pd.read_csv(paths["amr_calls"], dtype=str, keep_default_na=False)
    assert (calls["evidence_tier"] == "heuristic_motif").any()
    calls["evidence_tier"] = calls["evidence_tier"].str.upper()
    calls.to_csv(paths["amr_calls"], index=False)
    shouted = run_pipeline(
        **pipeline_kwargs(paths, output_dir=tmp_path / "u", config_overrides=CONFIG)
    )
    assert shouted["n_amr_calls_used"] == base["n_amr_calls_used"]


def test_out_of_range_dates_are_findings_not_crashes(tmp_path):
    samples = (
        "sample_id,organism,collection_date,country,source_sector\n"
        "A,E,1500-01-01,NO,human\nB,E,3000-01-01,NO,human\n"
    )
    ast = "sample_id,antimicrobial,result,method,interpretation_date\nA,cip,R,BMD,2999-01-01\n"
    _, _, codes = run(tmp_path, samples=samples, ast=ast)
    assert codes  # reported, no OutOfBoundsDatetime


def test_sub_annual_consumption_rows_are_not_duplicates(tmp_path):
    head = "year,country,care_sector,antimicrobial,ddd_total,population_n,period_days,source_name\n"
    two_months = (
        head + "2024,NO,community,J01,100,1000,31,SRC\n2024,NO,community,J01,90,1000,28,SRC\n"
    )
    _, _, codes = run(tmp_path, consumption=two_months)
    assert "duplicate_consumption_key" not in codes
    same_month = (
        head + "2024,NO,community,J01,100,1000,31,SRC\n2024,NO,community,J01,90,1000,31,SRC\n"
    )
    second = tmp_path / "second"
    second.mkdir()
    _, _, codes = run(second, consumption=same_month)
    assert "duplicate_consumption_key" in codes


def test_report_states_the_configured_confidence_level(tmp_path):
    ast = "sample_id,antimicrobial,result,method\nA,cip,R,BMD\n"
    out, _, _ = run(tmp_path, ast=ast, overrides={"analysis": {"confidence_level": 0.9}})
    text = (out / "report.md").read_text(encoding="utf-8")
    assert "Resistance (90% CI)" in text and "95% CI" not in text


def test_rows_dropped_from_context_streams_are_logged(tmp_path):
    inc = (
        "year,country,organism,antimicrobial,resistant_infection_n,population_n,source_name,case_definition\n"
        "2024,NO,E,cip,5,0,SRC,def\n2024,NO,E,cip,5,1000,SRC,def2\n"
    )
    out, _, _ = run(tmp_path, incidence=inc)
    log = pd.read_csv(out / "exclusion_log.csv")
    assert (log["reason"] == "invalid_row_excluded_from_context").sum() == 1


def test_feature_table_with_blank_feature_type_does_not_crash(tmp_path):
    (tmp_path / "f.csv").write_text(
        "gene,class,mechanism,feature_type,confidence\n"
        "qnrS1,QUINOLONE,ALLELEX,,1.0\nIncX,,,,0.9\ntetA,TET,ALLELEX,amr_gene,1.0\n",
        encoding="utf-8",
    )
    info = normalize_fishamr_feature_table(
        tmp_path / "f.csv", tmp_path / "c.csv", sample_id="ISO1", organism="E"
    )
    assert info["n_amr_calls"] == 1 and info["excluded_feature_types"] == {"missing": 2}


def test_infinite_and_implausible_values_are_errors(tmp_path):
    inc = (
        "year,country,organism,antimicrobial,resistant_infection_n,population_n,population_multiplier,source_name,case_definition\n"
        "2024,NO,E,cip,5,1000,inf,SRC,d1\n"
        "2024,NO,E,amp,5000,1000,100000,SRC,d2\n"
        "20244,NO,E,gen,5,1000,100000,SRC,d3\n"
    )
    cons = (
        "year,country,care_sector,antimicrobial,ddd_total,population_n,period_days,source_name\n"
        "2024,NO,community,J01,inf,1000,31,SRC\n2024,NO,hospital,J01,5,1000,100000,SRC\n"
    )
    out, _, codes = run(tmp_path, incidence=inc, consumption=cons)
    assert {
        "invalid_population_multiplier",
        "cases_exceed_population",
        "implausible_year",
        "non_finite_ddd_total",
        "implausible_period_days",
    } <= codes
    ctx = pd.read_csv(out / "incidence_context.csv")
    assert not ctx["incidence_per_population"].isin([float("inf")]).any()


def test_partial_interpretation_date_is_not_compared_with_collection_date(tmp_path):
    ast = (
        "sample_id,antimicrobial,result,method,interpretation_date\n"
        "A,cip,R,BMD,2024\nA,gen,R,BMD,2024-05\n"
    )
    _, _, codes = run(tmp_path, ast=ast)
    assert "interpretation_before_collection" not in codes
    assert "non_iso_interpretation_date" in codes


def test_zone_diameter_infinity_is_flagged(tmp_path):
    ast = "sample_id,antimicrobial,result,method,zone_mm\nA,cip,S,disk,inf\n"
    _, _, codes = run(tmp_path, ast=ast)
    assert "invalid_zone_mm" in codes


def test_binary_valid_accepts_one(tmp_path):
    (tmp_path / "t.csv").write_text(
        "isolate_id,antibiotic_name,probability_resistant,prediction_resistant,binary_valid\n"
        "A,cip,0.9,1,1\nB,cip,0.1,0,1\n",
        encoding="utf-8",
    )
    info = normalize_epiresnet_test_predictions(tmp_path / "t.csv", tmp_path / "o.csv")
    assert info["n_predictions_written"] == 2


def test_amrfinder_sample_id_keeps_dots_inside_the_file_name(tmp_path):
    report = tmp_path / "E.coli_1.tsv"
    report.write_text(
        "Gene symbol\tSequence name\tElement type\tMethod\tClass\tSubclass\n"
        "blaCTX-M-15\tbeta-lactamase\tAMR\tALLELEX\tBETA-LACTAM\tCEPHALOSPORIN\n",
        encoding="utf-8",
    )
    normalize_amrfinder_reports(
        [report],
        tmp_path / "c.csv",
        tmp_path / "s.csv",
        sample_ids=None,
        organism="E",
        tool_version="3.12",
        database_version="2024",
    )
    assert pd.read_csv(tmp_path / "c.csv")["sample_id"].iloc[0] == "E.coli_1"
