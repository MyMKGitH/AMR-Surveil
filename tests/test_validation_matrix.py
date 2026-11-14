"""One test per audit finding code: a valid base table, mutated just enough to trigger the code."""

import datetime as dt

import pandas as pd
import pytest

from amr_surveil.validation import audit_tables

REF = dt.date(2026, 11, 7)


def base():
    return {
        "samples": pd.DataFrame(
            {
                "sample_id": ["A", "B"],
                "organism": "E. coli",
                "collection_date": ["2024-01-01", "2024-02-01"],
                "country": "X",
                "source_sector": "human",
                "specimen": "blood",
                "sequencing_accession": "ACC",
                "source_study": "s",
                "data_license": "CC0",
            }
        ),
        "amr_calls": pd.DataFrame(
            {
                "sample_id": ["A"],
                "gene": "blaTEM-1",
                "identity_pct": "99",
                "coverage_pct": "100",
                "method": "AMRFinderPlus",
                "database": "D",
                "database_version": "1",
                "evidence_tier": "curated_database",
                "hit_method": "EXACTX",
            }
        ),
        "screening": pd.DataFrame(
            {
                "sample_id": ["A", "B"],
                "method": "AMRFinderPlus",
                "database": "D",
                "database_version": "1",
            }
        ),
        "ast": pd.DataFrame(
            {
                "sample_id": ["A"],
                "antimicrobial": "D",
                "result": "S",
                "mic": "1",
                "zone_mm": pd.NA,
                "method": "BMD",
                "breakpoint_standard": "EUCAST",
                "breakpoint_version": "14",
                "interpretation_date": "2025-01-01",
                "testing_lab": "L",
                "mic_unit": "mg/L",
            }
        ),
        "aggregate": pd.DataFrame(
            {
                "year": ["2024"],
                "country": "X",
                "organism": "E",
                "antimicrobial": "D",
                "resistant_n": ["1"],
                "tested_n": ["10"],
                "source_name": "s",
                "metric_definition": "m",
            }
        ),
        "predictions": pd.DataFrame(
            {
                "sample_id": ["A"],
                "antimicrobial": "D",
                "predicted_class": "predicted_resistant",
                "prediction_score": ["0.8"],
                "model_name": "M",
                "model_version": "1",
                "checkpoint_id": "c",
                "training_data_version": "t",
                "split_strategy": "s",
            }
        ),
        "incidence": pd.DataFrame(
            {
                "year": ["2024"],
                "country": "X",
                "organism": "E",
                "antimicrobial": "D",
                "resistant_infection_n": ["8"],
                "population_n": ["1000000"],
                "population_multiplier": ["100000"],
                "source_name": "s",
                "case_definition": "c",
            }
        ),
        "consumption": pd.DataFrame(
            {
                "year": ["2024"],
                "country": "X",
                "care_sector": "community",
                "antimicrobial": "Amoxicillin",
                "ddd_total": ["1000"],
                "population_n": ["100"],
                "period_days": ["365"],
                "source_name": "s",
                "data_source_type": "sales",
            }
        ),
    }


def audit(t):
    return audit_tables(
        t["samples"],
        t["amr_calls"],
        t["ast"],
        t["aggregate"],
        t["predictions"],
        pd.DataFrame(),
        t["incidence"],
        t["consumption"],
        screening=t["screening"],
        reference_date=REF,
    )


def test_base_tables_are_clean():
    out = audit(base())
    assert not set(out["severity"]) & {"ERROR", "WARNING"}, out[["code", "severity"]].to_string()


def with_(table, **cols):
    def mutate(t):
        for col, value in cols.items():
            t[table][col] = value
        return t

    return mutate


def dup_row(table, **overrides):
    def mutate(t):
        extra = t[table].iloc[[0]].assign(**overrides)
        t[table] = pd.concat([t[table], extra], ignore_index=True)
        return t

    return mutate


def replace(table, frame):
    def mutate(t):
        t[table] = frame
        return t

    return mutate


CASES = [
    ("empty_samples", "ERROR", replace("samples", pd.DataFrame(columns=["sample_id"]))),
    ("missing_sample_id", "ERROR", with_("samples", sample_id=[pd.NA, "B"])),
    ("nonstandard_source_sector", "WARNING", with_("samples", source_sector="spaceship")),
    ("missing_gene", "WARNING", with_("amr_calls", gene=pd.NA)),
    ("invalid_identity_pct", "ERROR", with_("amr_calls", identity_pct="150")),
    ("invalid_coverage_pct", "ERROR", with_("amr_calls", coverage_pct="-3")),
    ("duplicate_amr_call", "WARNING", dup_row("amr_calls")),
    ("evidence_tier_unspecified", "INFO", with_("amr_calls", evidence_tier=pd.NA)),
    (
        "calls_for_unscreened_sample",
        "WARNING",
        replace(
            "screening",
            pd.DataFrame(
                {"sample_id": ["B"], "method": "m", "database": "d", "database_version": "1"}
            ),
        ),
    ),
    ("multiple_screening_records", "WARNING", dup_row("screening", database_version="2")),
    ("invalid_zone_mm", "WARNING", with_("ast", zone_mm="-5")),
    ("conflicting_duplicate_ast", "ERROR", dup_row("ast", result="R")),
    ("duplicate_ast_record", "WARNING", dup_row("ast")),
    ("non_iso_interpretation_date", "WARNING", with_("ast", interpretation_date="05/06/2025")),
    ("invalid_aggregate_year", "WARNING", with_("aggregate", year="20x4")),
    ("invalid_aggregate_counts", "ERROR", with_("aggregate", resistant_n="11")),
    ("non_integer_aggregate_counts", "ERROR", with_("aggregate", resistant_n="1.5")),
    ("zero_aggregate_denominator", "WARNING", with_("aggregate", resistant_n="0", tested_n="0")),
    ("missing_source_name", "WARNING", with_("aggregate", source_name=pd.NA)),
    (
        "missing_prediction_value",
        "WARNING",
        with_("predictions", predicted_mic=pd.NA, predicted_class=pd.NA),
    ),
    ("non_probability_prediction_score", "WARNING", with_("predictions", prediction_score="1.5")),
    ("zero_population_denominator", "ERROR", with_("incidence", population_n="0")),
    ("negative_resistant_infection_n", "ERROR", with_("incidence", resistant_infection_n="-1")),
    ("non_integer_resistant_infection_n", "ERROR", with_("incidence", resistant_infection_n="2.5")),
    ("invalid_population_multiplier", "ERROR", with_("incidence", population_multiplier="0")),
    ("zero_period_days", "ERROR", with_("consumption", period_days="0")),
    ("negative_ddd_total", "ERROR", with_("consumption", ddd_total="-5")),
    ("nonstandard_care_sector", "WARNING", with_("consumption", care_sector="spa")),
    (
        "nonstandard_consumption_source_type",
        "WARNING",
        with_("consumption", data_source_type="guess"),
    ),
]


@pytest.mark.parametrize("code,severity,mutate", CASES, ids=[c[0] for c in CASES])
def test_finding_is_raised_with_expected_severity(code, severity, mutate):
    out = audit(mutate(base()))
    hit = out[out.code == code]
    assert not hit.empty, f"{code} not raised; got {sorted(out.code)}"
    assert set(hit.severity) == {severity}
    if code != "empty_samples":  # table-level finding has no row count
        assert (hit.row_count >= 1).all()
