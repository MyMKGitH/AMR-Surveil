import pytest

from amr_surveil.config import DEFAULTS, load_config, resolve_reference_date, validate_config
from amr_surveil.exceptions import ConfigurationError


def test_defaults_are_valid_and_conservative():
    cfg = load_config(None)
    assert cfg["quality"]["strict"] is True
    assert set(cfg["analysis"]["excluded_evidence_tiers"]) == {"heuristic_motif", "homology_search"}
    assert cfg["study"]["predictions_independence_attested"] is False
    assert cfg["study"]["sample_selection_documented"] is False


@pytest.mark.parametrize(
    "loaded,match",
    [
        ({"bogus": {}}, "Unsupported configuration section"),
        ({"analysis": {"typo": 1}}, "Unsupported analysis setting"),
        ({"analysis": {"confidence_level": 1.5}}, "between 0 and 1"),
        ({"analysis": {"confidence_level": True}}, "number"),
        ({"analysis": {"small_denominator_warning_threshold": 2.5}}, "whole number"),
        ({"analysis": {"small_denominator_warning_threshold": 0}}, "whole number"),
        ({"analysis": {"deduplication_rule": "magic"}}, "Unsupported deduplication_rule"),
        ({"analysis": {"excluded_evidence_tiers": ["nope"]}}, "Unknown evidence tier"),
        ({"analysis": {"reference_date": "07/11/2026"}}, "ISO date"),
        ({"quality": {"strict": "false"}}, "YAML boolean"),
        ({"study": {"sample_selection_documented": "yes"}}, "YAML boolean"),
        ({"study": {"name": 5}}, "must be text"),
        ({"analysis": "x"}, "must be an object"),
    ],
)
def test_invalid_config_rejected(loaded, match):
    with pytest.raises(ConfigurationError, match=match):
        validate_config(loaded)


def test_reference_date_resolution_and_no_default_mutation():
    cfg = validate_config({"analysis": {"reference_date": "2026-01-31"}})
    assert resolve_reference_date(cfg).isoformat() == "2026-01-31"
    assert resolve_reference_date(load_config(None)) is not None
    cfg["analysis"]["excluded_evidence_tiers"].append("x")
    assert "x" not in DEFAULTS["analysis"]["excluded_evidence_tiers"]


def test_missing_config_file(tmp_path):
    with pytest.raises(ConfigurationError, match="Could not read"):
        load_config(tmp_path / "nope.yaml")
