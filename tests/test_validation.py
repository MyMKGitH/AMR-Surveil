import datetime as dt

import pandas as pd
import pytest

from amr_surveil.constants import ISSUE_COLUMNS
from amr_surveil.validation import audit_tables

REF = dt.date(2026, 11, 7)
E = pd.DataFrame()


def audit(samples=None, calls=None, ast=None, screening=None, **kw):
    samples = (
        samples
        if samples is not None
        else pd.DataFrame(
            {
                "sample_id": ["A", "B"],
                "organism": ["E. coli"] * 2,
                "collection_date": ["2024-01-01", "2024-02-01"],
            }
        )
    )
    return audit_tables(
        samples,
        calls if calls is not None else E,
        ast if ast is not None else E,
        E,
        E,
        E,
        screening=screening,
        reference_date=REF,
        **kw,
    )


def codes(frame, severity=None):
    f = frame if severity is None else frame[frame.severity == severity]
    return set(f["code"])


def test_clean_minimal_input_has_no_errors():
    out = audit()
    assert list(out.columns) == ISSUE_COLUMNS
    assert "ERROR" not in set(out["severity"])


@pytest.mark.parametrize(
    "date,code,severity",
    [
        ("03/05/2024", "non_iso_collection_date", "WARNING"),
        ("2026-11-08", "future_collection_date", "ERROR"),  # one day after the reference date
        ("2024-02-30", "invalid_collection_date", None),
    ],
)
def test_collection_date_rules(date, code, severity):
    s = pd.DataFrame({"sample_id": ["A"], "collection_date": [date]})
    out = audit(samples=s)
    assert code in codes(out)
    if severity:
        assert code in codes(out, severity)


def test_reference_date_is_inclusive():
    s = pd.DataFrame({"sample_id": ["A"], "collection_date": ["2026-11-07"]})
    assert "future_collection_date" not in codes(audit(samples=s))


def test_duplicate_and_label_variant_findings():
    s = pd.DataFrame({"sample_id": ["A", "A"], "organism": ["E. coli", "e. coli "]})
    out = audit(samples=s)
    assert "duplicate_sample_id" in codes(out, "ERROR")
    assert "inconsistent_organism_labels" in codes(out, "WARNING")


def ast_frame(**over):
    base = {
        "sample_id": "A",
        "antimicrobial": "D",
        "result": "S",
        "mic": "1",
        "method": "BMD",
        "breakpoint_standard": "EUCAST",
        "breakpoint_version": "14",
        "interpretation_date": "2025-01-01",
    }
    return pd.DataFrame([{**base, **over}])


def test_blank_result_is_warning_but_wrong_category_is_error():
    assert "missing_result" in codes(audit(ast=ast_frame(result=pd.NA)), "WARNING")
    assert "invalid_ast_category" not in codes(audit(ast=ast_frame(result=pd.NA)))
    assert "invalid_ast_category" in codes(audit(ast=ast_frame(result="Intermediate")), "ERROR")


def test_mic_findings():
    assert "invalid_mic" in codes(audit(ast=ast_frame(mic="1,5")))
    assert "nonpositive_mic" in codes(audit(ast=ast_frame(mic="0")))
    assert "censored_mic_present" in codes(audit(ast=ast_frame(mic="<=0.5")), "INFO")
    assert "invalid_mic" not in codes(audit(ast=ast_frame(mic=">=32")))


def test_unmatched_ast_sample_and_interpretation_before_collection():
    assert "unmatched_sample_id" in codes(audit(ast=ast_frame(sample_id="ZZZ")))
    assert "interpretation_before_collection" in codes(
        audit(ast=ast_frame(interpretation_date="2023-01-01"))
    )


def call_frame(**over):
    base = {
        "sample_id": "A",
        "gene": "blaTEM-1",
        "identity_pct": "99",
        "coverage_pct": "100",
        "evidence_tier": "curated_database",
        "hit_method": "EXACTX",
        "method": "AMRFinderPlus",
        "database": "D",
        "database_version": "1",
    }
    return pd.DataFrame([{**base, **over}])


def test_amr_call_findings():
    assert "invalid_evidence_tier" in codes(audit(calls=call_frame(evidence_tier="magic")), "ERROR")
    assert "low_confidence_evidence_tier" in codes(
        audit(calls=call_frame(evidence_tier="heuristic_motif"))
    )
    assert "low_identity_or_coverage_review" in codes(audit(calls=call_frame(identity_pct="70")))
    assert "low_identity_or_coverage_review" in codes(audit(calls=call_frame(coverage_pct="40")))
    assert "partial_or_truncated_hit" in codes(
        audit(calls=call_frame(hit_method="PARTIALX"))
    ) or "partial_or_internal_stop_hit" in codes(audit(calls=call_frame(hit_method="PARTIALX")))


def test_screening_provenance_findings():
    scr = pd.DataFrame(
        {
            "sample_id": ["A"],
            "method": ["AMRFinderPlus"],
            "database": [pd.NA],
            "database_version": [pd.NA],
        }
    )
    out = audit(screening=scr)
    assert any(c.startswith("missing_screening_") for c in codes(out))


def test_duplicate_aggregate_key_is_error():
    agg = pd.DataFrame(
        {
            "year": ["2024"] * 2,
            "country": ["X"] * 2,
            "organism": ["E"] * 2,
            "antimicrobial": ["D"] * 2,
            "resistant_n": ["1", "2"],
            "tested_n": ["10", "10"],
            "source_name": ["s"] * 2,
            "metric_definition": ["m"] * 2,
        }
    )
    out = audit_tables(pd.DataFrame({"sample_id": ["A"]}), E, E, agg, E, E, reference_date=REF)
    assert "duplicate_aggregate_key" in codes(out, "ERROR")


def test_findings_sorted_by_severity():
    out = audit(
        samples=pd.DataFrame({"sample_id": ["A", "A"], "collection_date": ["x", "2024-01-01"]})
    )
    order = {"ERROR": 0, "WARNING": 1, "INFO": 2}
    ranks = [order[s] for s in out["severity"]]
    assert ranks == sorted(ranks)
