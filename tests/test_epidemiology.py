import pandas as pd
import pytest

from amr_surveil.epidemiology import resistance_indicators, wilson_interval


def test_wilson_edges_and_invalid_inputs():
    low, high = wilson_interval(0, 10)
    assert low == pytest.approx(0.0)
    assert high > 0
    low, high = wilson_interval(10, 10)
    assert high == pytest.approx(100.0)
    assert low < 100
    with pytest.raises(ValueError):
        wilson_interval(11, 10)
    with pytest.raises(ValueError):
        wilson_interval(1, 0)


def test_resistance_denominator_keeps_i_but_r_is_numerator():
    samples = pd.DataFrame(
        {
            "sample_id": ["A", "B", "C"],
            "organism": ["E. coli"] * 3,
            "collection_date": ["2024-01-01", "2024-01-02", "2024-01-03"],
            "country": ["X"] * 3,
            "source_sector": ["human"] * 3,
        }
    )
    ast = pd.DataFrame(
        {
            "sample_id": ["A", "B", "C"],
            "organism": ["E. coli"] * 3,
            "antimicrobial": ["DrugA"] * 3,
            "result": ["S", "I", "R"],
            "method": ["BMD"] * 3,
            "breakpoint_standard": ["EUCAST"] * 3,
            "breakpoint_version": ["16.1"] * 3,
        }
    )
    result = resistance_indicators(samples, ast)
    row = result.iloc[0]
    assert row.n_tested == 3
    assert row.n_resistant == 1
    assert row.resistant_pct == pytest.approx(100 / 3, abs=0.001)


def test_breakpoint_context_is_retained_in_groups():
    samples = pd.DataFrame(
        {
            "sample_id": ["A", "B"],
            "organism": ["E. coli"] * 2,
            "collection_date": ["2024-01-01", "2024-01-02"],
            "country": ["X"] * 2,
            "source_sector": ["human"] * 2,
        }
    )
    ast = pd.DataFrame(
        {
            "sample_id": ["A", "B"],
            "organism": ["E. coli"] * 2,
            "antimicrobial": ["DrugA"] * 2,
            "result": ["R", "S"],
            "method": ["BMD", "Disk"],
            "breakpoint_standard": ["EUCAST"] * 2,
            "breakpoint_version": ["16.1"] * 2,
        }
    )
    result = resistance_indicators(samples, ast)
    assert len(result) == 2
    assert set(result["method"]) == {"BMD", "Disk"}


def test_duplicate_manifest_keys_are_excluded_in_exploratory_indicator_run():
    samples = pd.DataFrame(
        {
            "sample_id": ["A", "A"],
            "organism": ["E. coli", "E. coli"],
            "collection_date": ["2024-01-01", "2024-01-01"],
            "country": ["X", "X"],
            "source_sector": ["human", "human"],
        }
    )
    ast = pd.DataFrame({"sample_id": ["A"], "antimicrobial": ["DrugA"], "result": ["R"]})
    assert resistance_indicators(samples, ast).empty


def test_patient_deduplication_requires_real_patient_ids():
    from amr_surveil.exceptions import ConfigurationError

    samples = pd.DataFrame(
        {
            "sample_id": ["A"],
            "patient_id": [pd.NA],
            "organism": ["E. coli"],
            "collection_date": ["2024-01-01"],
        }
    )
    ast = pd.DataFrame({"sample_id": ["A"], "result": ["R"]})
    with pytest.raises(ConfigurationError, match="at least one populated patient_id"):
        resistance_indicators(samples, ast, deduplication_rule="first_per_patient_species_year")


def test_patient_dedup_preserves_rows_without_patient_identifier():
    samples = pd.DataFrame(
        {
            "sample_id": ["A", "B", "C"],
            "patient_id": ["P1", "P1", pd.NA],
            "organism": ["E. coli"] * 3,
            "collection_date": ["2024-01-01", "2024-02-01", "2024-03-01"],
        }
    )
    ast = pd.DataFrame(
        {"sample_id": ["A", "B", "C"], "antimicrobial": ["DrugA"] * 3, "result": ["R"] * 3}
    )
    result = resistance_indicators(
        samples, ast, deduplication_rule="first_per_patient_species_year"
    )
    # A is the first P1 isolate; C is not a patient-linked record and is retained.
    assert result["n_tested"].sum() == 2


# ---- reference values: Newcombe (1998) Stat Med 17:857-872, Table I, Wilson score method (95%)
@pytest.mark.parametrize(
    "x,n,low,high",
    [(81, 263, 25.53, 36.62), (0, 20, 0.0, 16.11), (20, 20, 83.89, 100.0)],
)
def test_wilson_matches_published_reference_values(x, n, low, high):
    got_low, got_high = wilson_interval(x, n)
    assert got_low == pytest.approx(low, abs=0.01)
    assert got_high == pytest.approx(high, abs=0.01)


def test_wilson_independent_closed_form_and_confidence_level():
    z95, z90 = 1.959964, 1.644854
    for x, n, z in ((8, 10, z95), (50, 100, z95), (3, 40, z90)):
        p = x / n
        centre = (p + z * z / (2 * n)) / (1 + z * z / n)
        half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / (1 + z * z / n)
        level = 0.95 if z == z95 else 0.90
        low, high = wilson_interval(x, n, level)
        assert low == pytest.approx(100 * (centre - half), abs=1e-3)
        assert high == pytest.approx(100 * (centre + half), abs=1e-3)
    assert (
        wilson_interval(8, 10, 0.90)[0] > wilson_interval(8, 10, 0.95)[0]
    )  # narrower at lower confidence


def test_small_denominator_flag_boundary_is_strict_less_than():
    def table(n):
        samples = pd.DataFrame(
            {
                "sample_id": [f"S{i}" for i in range(n)],
                "organism": "E",
                "collection_date": "2024-01-01",
            }
        )
        ast = pd.DataFrame({"sample_id": samples.sample_id, "antimicrobial": "D", "result": "S"})
        return resistance_indicators(samples, ast, low_denominator_threshold=10)

    assert bool(table(9).small_denominator_warning.iloc[0]) is True
    assert (
        bool(table(10).small_denominator_warning.iloc[0]) is False
    )  # n == threshold is not "small"
    assert bool(table(11).small_denominator_warning.iloc[0]) is False
