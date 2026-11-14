"""Population-based incidence and antimicrobial-consumption indicators."""

from __future__ import annotations

import numpy as np
import pandas as pd


def population_incidence_context(incidence: pd.DataFrame) -> pd.DataFrame:
    """Calculate cases per stated population multiplier from source-defined counts.

    The caller must provide the case definition and matching population denominator.
    Counts are never inferred from isolate-level AST or aggregate resistance proportions.
    """
    if incidence.empty:
        return pd.DataFrame(
            columns=list(incidence.columns) + ["incidence_per_population", "evidence_level"]
        )
    frame = incidence.copy()
    frame["resistant_infection_n"] = pd.to_numeric(frame["resistant_infection_n"], errors="coerce")
    frame["population_n"] = pd.to_numeric(frame["population_n"], errors="coerce")
    if "population_multiplier" not in frame.columns:
        frame["population_multiplier"] = 100_000
    else:
        raw_multiplier = frame["population_multiplier"]
        parsed_multiplier = pd.to_numeric(raw_multiplier, errors="coerce")
        frame["population_multiplier"] = parsed_multiplier.where(raw_multiplier.notna(), 100_000)
    valid = (
        frame["resistant_infection_n"].notna()
        & frame["population_n"].notna()
        & frame["population_multiplier"].notna()
        & (frame["resistant_infection_n"] >= 0)
        & (frame["population_n"] > 0)
        & (frame["population_multiplier"] > 0)
        & np.isfinite(
            frame[["resistant_infection_n", "population_n", "population_multiplier"]]
        ).all(axis=1)
    )
    frame = frame[valid].copy()
    frame["incidence_per_population"] = (
        frame["resistant_infection_n"] / frame["population_n"] * frame["population_multiplier"]
    ).round(6)
    frame["evidence_level"] = "population-based incidence; aggregate and not isolate-linked"
    return frame


def antimicrobial_consumption_context(consumption: pd.DataFrame) -> pd.DataFrame:
    """Calculate DDD/1,000 inhabitants/day from supplied DDD, population and days.

    This is a standardised consumption indicator, not a patient-level exposure or
    ingestion measure. Input rows are retained separately; no implicit summation occurs.
    """
    if consumption.empty:
        return pd.DataFrame(
            columns=list(consumption.columns)
            + ["ddd_per_1000_inhabitants_per_day", "evidence_level"]
        )
    frame = consumption.copy()
    for col in ("ddd_total", "population_n", "period_days"):
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    valid = (
        frame["ddd_total"].notna()
        & frame["population_n"].notna()
        & frame["period_days"].notna()
        & (frame["ddd_total"] >= 0)
        & (frame["population_n"] > 0)
        & (frame["period_days"] > 0)
        & np.isfinite(frame[["ddd_total", "population_n", "period_days"]]).all(axis=1)
    )
    frame = frame[valid].copy()
    frame["ddd_per_1000_inhabitants_per_day"] = (
        frame["ddd_total"] / (frame["population_n"] * frame["period_days"]) * 1000
    ).round(6)
    frame["evidence_level"] = "aggregate antimicrobial consumption; not patient-level use"
    return frame
