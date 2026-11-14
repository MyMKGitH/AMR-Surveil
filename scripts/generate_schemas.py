"""Write schemas/*.schema.json from the column contracts in amr_surveil.constants.

Run ``python scripts/generate_schemas.py`` after changing a contract; ``tests/test_schema_contracts.py``
fails when the checked-in schemas are out of date.
"""

from __future__ import annotations

import json
from pathlib import Path

from amr_surveil import constants as C

BASE = "https://github.com/MyMKGitH/amr-surveil/schemas"
STR = {"type": ["string", "null"]}
SPECS: dict[str, dict] = {
    "samples": {
        "title": "Sample / isolate manifest row",
        "columns": C.SAMPLE_COLUMNS,
        "required": ["sample_id"],
        "props": {
            "sample_id": {"type": "string", "minLength": 1},
            "collection_date": {
                "type": ["string", "null"],
                "description": "ISO 8601: YYYY, YYYY-MM or YYYY-MM-DD",
            },
            "source_sector": {"type": ["string", "null"], "examples": sorted(C.VALID_SECTORS)},
        },
    },
    "amr_calls": {
        "title": "AMR determinant call (positive hits only)",
        "columns": C.AMR_CALL_COLUMNS,
        "required": ["sample_id", "gene"],
        "props": {
            "sample_id": {"type": "string", "minLength": 1},
            "gene": {"type": "string", "minLength": 1},
            "identity_pct": {"type": ["number", "string", "null"], "description": "0-100"},
            "coverage_pct": {"type": ["number", "string", "null"], "description": "0-100"},
            "evidence_tier": {"type": ["string", "null"], "enum": [*C.EVIDENCE_TIERS, None]},
        },
    },
    "screening": {
        "title": "Screening provenance: which samples were genotyped, with what",
        "columns": C.SCREENING_COLUMNS,
        "required": ["sample_id"],
        "props": {"sample_id": {"type": "string", "minLength": 1}},
    },
    "ast": {
        "title": "Phenotypic AST record",
        "columns": C.AST_COLUMNS,
        "required": ["sample_id", "antimicrobial", "result"],
        "props": {
            "sample_id": {"type": "string", "minLength": 1},
            "antimicrobial": {"type": "string", "minLength": 1},
            "result": {"type": ["string", "null"], "enum": ["S", "I", "R", None]},
            "mic": {
                "type": ["number", "string", "null"],
                "description": "Optional operator <=,<,>=,> then a number; no decimal commas",
            },
            "zone_mm": {"type": ["number", "string", "null"]},
        },
    },
    "breakpoints": {
        "title": "Breakpoint row used to re-check recorded categories",
        "columns": C.BREAKPOINT_COLUMNS,
        "required": ["organism", "antimicrobial", "standard", "version"],
        "props": {},
    },
    "aggregate_surveillance": {
        "title": "Aggregate resistance counts (not isolate-linked)",
        "columns": C.AGGREGATE_COLUMNS,
        "required": [
            "year",
            "country",
            "organism",
            "antimicrobial",
            "resistant_n",
            "tested_n",
            "source_name",
            "metric_definition",
        ],
        "props": {},
    },
    "incidence": {
        "title": "Population-based resistant-infection counts",
        "columns": C.INCIDENCE_COLUMNS,
        "required": [
            "year",
            "country",
            "organism",
            "antimicrobial",
            "resistant_infection_n",
            "population_n",
            "source_name",
            "case_definition",
        ],
        "props": {},
    },
    "consumption": {
        "title": "Antimicrobial consumption (DDD)",
        "columns": C.CONSUMPTION_COLUMNS,
        "required": [
            "year",
            "country",
            "care_sector",
            "antimicrobial",
            "ddd_total",
            "population_n",
            "period_days",
            "source_name",
        ],
        "props": {},
    },
    "model_predictions": {
        "title": "Model prediction (never observed AST)",
        "columns": C.PREDICTION_COLUMNS,
        "required": ["sample_id", "antimicrobial"],
        "props": {"sample_id": {"type": "string", "minLength": 1}},
    },
    "plasmid_calls": {
        "title": "Plasmid / replicon evidence",
        "columns": C.PLASMID_COLUMNS,
        "required": ["sample_id", "replicon"],
        "props": {"evidence_tier": {"type": ["string", "null"], "enum": [*C.EVIDENCE_TIERS, None]}},
    },
    "gene_drug_map": {
        "title": "Cited gene-to-antimicrobial mapping",
        "columns": C.GENE_MAP_COLUMNS,
        "required": ["organism", "gene", "antimicrobial", "mapping_source"],
        "props": {"match_type": {"type": ["string", "null"], "enum": ["exact", "prefix", None]}},
    },
}


def build(name: str, spec: dict) -> dict:
    required_default = {"type": "string", "minLength": 1}
    props = {}
    for column in spec["columns"]:
        default = required_default if column in spec["required"] else STR
        props[column] = spec["props"].get(column, default)
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": f"{BASE}/{name}.schema.json",
        "title": f"AMR-Surveil: {spec['title']}",
        "type": "object",
        "required": spec["required"],
        "properties": props,
        "additionalProperties": True,
    }


def main() -> int:
    out = Path(__file__).resolve().parents[1] / "schemas"
    out.mkdir(exist_ok=True)
    for name, spec in SPECS.items():
        (out / f"{name}.schema.json").write_text(
            json.dumps(build(name, spec), indent=2) + "\n", encoding="utf-8"
        )
        print("wrote", out / f"{name}.schema.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
