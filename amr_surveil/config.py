"""Validated run configuration.

Every setting that can change a result is declared here, validated, echoed into the run
manifest and documented in ``docs/METHODS.md``. Unknown sections or keys are rejected so a
typo can never silently fall back to a default.
"""

from __future__ import annotations

import copy
import datetime as dt
from pathlib import Path
from typing import Any

import yaml

from amr_surveil.constants import (
    DEDUPLICATION_RULES,
    DEFAULT_EXCLUDED_EVIDENCE_TIERS,
    EVIDENCE_TIERS,
)
from amr_surveil.exceptions import ConfigurationError

DEFAULTS: dict[str, dict[str, Any]] = {
    "analysis": {
        "deduplication_rule": "none",
        "confidence_level": 0.95,
        "small_denominator_warning_threshold": 10,
        # ``null`` => today's UTC date, which is recorded in the run manifest.
        "reference_date": None,
        "excluded_evidence_tiers": list(DEFAULT_EXCLUDED_EVIDENCE_TIERS),
        "review_identity_pct": 90.0,
        "review_coverage_pct": 60.0,
        "concordance_min_isolates": 10,
        "testing_volume_jump_ratio": 2.0,
    },
    "quality": {"strict": True},
    "study": {
        "name": "",
        "question": "",
        "inclusion_criteria": "",
        "exclusion_criteria": "",
        "sample_selection_documented": False,
        "selection_description": "",
        "predictions_independence_attested": False,
    },
}


def _require_bool(section: str, key: str, value: Any) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{section}.{key} must be a YAML boolean (true or false)")
    return value


def _require_text(section: str, key: str, value: Any) -> str:
    if value is None:
        return ""
    if not isinstance(value, str):
        raise ConfigurationError(f"{section}.{key} must be text")
    return value.strip()


def _require_whole(section: str, key: str, value: Any, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or float(value) != int(value):
        raise ConfigurationError(f"{section}.{key} must be a whole number >= {minimum}")
    if int(value) < minimum:
        raise ConfigurationError(f"{section}.{key} must be a whole number >= {minimum}")
    return int(value)


def _require_number(section: str, key: str, value: Any, low: float, high: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{section}.{key} must be a number")
    number = float(value)
    if not low <= number <= high:
        raise ConfigurationError(f"{section}.{key} must be between {low} and {high}")
    return number


def validate_config(loaded: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Merge ``loaded`` over the defaults and validate every value."""
    if not isinstance(loaded, dict):
        raise ConfigurationError("Configuration root must be a YAML object")
    unknown_sections = set(loaded) - set(DEFAULTS)
    if unknown_sections:
        raise ConfigurationError(
            f"Unsupported configuration section(s): {', '.join(sorted(unknown_sections))}"
        )
    config = copy.deepcopy(DEFAULTS)
    for section, defaults in DEFAULTS.items():
        value = loaded.get(section) or {}
        if not isinstance(value, dict):
            raise ConfigurationError(f"Configuration section '{section}' must be an object")
        unknown = set(value) - set(defaults)
        if unknown:
            raise ConfigurationError(
                f"Unsupported {section} setting(s): {', '.join(sorted(unknown))}"
            )
        config[section].update(value)

    a = config["analysis"]
    if a["deduplication_rule"] not in DEDUPLICATION_RULES:
        raise ConfigurationError(
            f"Unsupported deduplication_rule: {a['deduplication_rule']} "
            f"(choose from {', '.join(DEDUPLICATION_RULES)})"
        )
    confidence = a["confidence_level"]
    if isinstance(confidence, bool) or not isinstance(confidence, (int, float)):
        raise ConfigurationError("analysis.confidence_level must be a number between 0 and 1")
    if not 0 < float(confidence) < 1:
        raise ConfigurationError("analysis.confidence_level must be between 0 and 1")
    a["confidence_level"] = float(confidence)
    a["small_denominator_warning_threshold"] = _require_whole(
        "analysis",
        "small_denominator_warning_threshold",
        a["small_denominator_warning_threshold"],
        1,
    )
    a["concordance_min_isolates"] = _require_whole(
        "analysis", "concordance_min_isolates", a["concordance_min_isolates"], 1
    )
    a["review_identity_pct"] = _require_number(
        "analysis", "review_identity_pct", a["review_identity_pct"], 0, 100
    )
    a["review_coverage_pct"] = _require_number(
        "analysis", "review_coverage_pct", a["review_coverage_pct"], 0, 100
    )
    ratio = _require_number(
        "analysis", "testing_volume_jump_ratio", a["testing_volume_jump_ratio"], 1.0, 1e6
    )
    a["testing_volume_jump_ratio"] = ratio

    tiers = a["excluded_evidence_tiers"]
    if not isinstance(tiers, list) or any(not isinstance(t, str) for t in tiers):
        raise ConfigurationError("analysis.excluded_evidence_tiers must be a list of tier names")
    bad = sorted(set(tiers) - set(EVIDENCE_TIERS))
    if bad:
        raise ConfigurationError(
            f"Unknown evidence tier(s) in analysis.excluded_evidence_tiers: {', '.join(bad)}"
        )
    a["excluded_evidence_tiers"] = sorted(set(tiers))

    ref = a["reference_date"]
    if ref is not None:
        try:
            a["reference_date"] = dt.date.fromisoformat(str(ref)).isoformat()
        except ValueError as exc:
            raise ConfigurationError(
                "analysis.reference_date must be an ISO date (YYYY-MM-DD) or null"
            ) from exc

    config["quality"]["strict"] = _require_bool("quality", "strict", config["quality"]["strict"])

    s = config["study"]
    for key in (
        "name",
        "question",
        "inclusion_criteria",
        "exclusion_criteria",
        "selection_description",
    ):
        s[key] = _require_text("study", key, s[key])
    for key in ("sample_selection_documented", "predictions_independence_attested"):
        s[key] = _require_bool("study", key, s[key])
    return config


def load_config(path: str | Path | None) -> dict[str, dict[str, Any]]:
    """Load a YAML configuration (or return validated defaults when ``path`` is None)."""
    if path is None:
        return validate_config({})
    try:
        loaded = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigurationError(f"Could not read configuration {path}: {exc}") from exc
    return validate_config(loaded)


def resolve_reference_date(config: dict[str, dict[str, Any]]) -> dt.date:
    """The date against which 'future' collection dates are judged (recorded in the manifest)."""
    ref = config["analysis"]["reference_date"]
    if ref:
        return dt.date.fromisoformat(ref)
    return dt.datetime.now(dt.timezone.utc).date()
