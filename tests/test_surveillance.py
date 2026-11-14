import pandas as pd
import pytest

from amr_surveil.surveillance import (
    antimicrobial_consumption_context,
    population_incidence_context,
)


def test_population_incidence_uses_explicit_denominator_and_multiplier():
    frame = pd.DataFrame(
        {
            "resistant_infection_n": [8],
            "population_n": [1_000_000],
            "population_multiplier": [100_000],
            "case_definition": ["resistant bloodstream infection episodes"],
        }
    )
    result = population_incidence_context(frame)
    assert result.loc[0, "incidence_per_population"] == pytest.approx(0.8)
    assert result.loc[0, "evidence_level"].startswith("population-based incidence")


def test_population_incidence_default_multiplier_is_explicit_in_output():
    frame = pd.DataFrame({"resistant_infection_n": [5], "population_n": [100_000]})
    result = population_incidence_context(frame)
    assert result.loc[0, "population_multiplier"] == 100_000
    assert result.loc[0, "incidence_per_population"] == pytest.approx(5.0)


def test_ddd_consumption_formula_uses_population_and_period():
    frame = pd.DataFrame(
        {"ddd_total": [1_000_000], "population_n": [100_000], "period_days": [365]}
    )
    result = antimicrobial_consumption_context(frame)
    assert result.loc[0, "ddd_per_1000_inhabitants_per_day"] == pytest.approx(27.39726, abs=1e-5)
    assert "not patient-level" in result.loc[0, "evidence_level"]


def test_invalid_rows_are_not_calculated():
    incidence = pd.DataFrame({"resistant_infection_n": [2], "population_n": [0]})
    consumption = pd.DataFrame({"ddd_total": [1], "population_n": [100], "period_days": [0]})
    assert population_incidence_context(incidence).empty
    assert antimicrobial_consumption_context(consumption).empty
