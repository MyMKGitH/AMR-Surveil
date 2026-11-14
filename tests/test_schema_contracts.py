import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

jsonschema = pytest.importorskip("jsonschema")
ROOT = Path(__file__).resolve().parents[1]

DEMO_FILES = {
    "samples": "samples.csv",
    "amr_calls": "amr_calls.csv",
    "screening": "screening.csv",
    "ast": "ast.csv",
    "breakpoints": "breakpoints.csv",
    "gene_drug_map": "gene_drug_map.csv",
    "aggregate_surveillance": "aggregate_surveillance.csv",
    "model_predictions": "model_predictions.csv",
    "plasmid_calls": "plasmid_calls.csv",
    "incidence": "incidence.csv",
    "consumption": "consumption.csv",
}


def _generator():
    spec = importlib.util.spec_from_file_location(
        "generate_schemas", ROOT / "scripts" / "generate_schemas.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_checked_in_schemas_match_the_column_contracts():
    gen = _generator()
    for name, spec in gen.SPECS.items():
        on_disk = json.loads((ROOT / "schemas" / f"{name}.schema.json").read_text(encoding="utf-8"))
        assert on_disk == gen.build(name, spec), f"{name}: run python scripts/generate_schemas.py"


@pytest.mark.parametrize("table,filename", DEMO_FILES.items())
def test_every_demo_row_validates_against_its_schema(table, filename):
    schema = json.loads((ROOT / "schemas" / f"{table}.schema.json").read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator.check_schema(schema)
    frame = pd.read_csv(ROOT / "data" / "demo" / filename, dtype=str, keep_default_na=False)
    assert set(frame.columns) <= set(schema["properties"]), table
    assert set(schema["required"]) <= set(frame.columns), table
    validator = jsonschema.Draft202012Validator(schema)
    for row in frame.to_dict("records"):
        record = {k: (None if v == "" else v) for k, v in row.items()}
        errors = list(validator.iter_errors(record))
        assert not errors, (table, record, errors[0].message)


def test_schema_rejects_a_bad_ast_category():
    schema = json.loads((ROOT / "schemas" / "ast.schema.json").read_text(encoding="utf-8"))
    assert list(
        jsonschema.Draft202012Validator(schema).iter_errors(
            {"sample_id": "A", "antimicrobial": "D", "result": "X"}
        )
    )


def test_example_adapter_mappings_are_valid_json_with_supported_targets():
    from amr_surveil.constants import AMR_CALL_COLUMNS, PLASMID_COLUMNS, PREDICTION_COLUMNS

    allowed = {
        "fishamr_mapping.example.json": AMR_CALL_COLUMNS,
        "fishamr_plasmid_mapping.example.json": PLASMID_COLUMNS,
        "epiresnet_mapping.example.json": PREDICTION_COLUMNS,
    }
    for name, columns in allowed.items():
        mapping = json.loads((ROOT / "config" / name).read_text(encoding="utf-8"))
        assert set(mapping["columns"]) <= set(columns), name
