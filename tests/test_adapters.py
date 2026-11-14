import json
from pathlib import Path

import pandas as pd
import pytest

from amr_surveil.adapters.epiresnet import normalize_epiresnet_predictions
from amr_surveil.adapters.fishamr import normalize_fishamr_calls
from amr_surveil.exceptions import ConfigurationError


def test_fishamr_adapter_maps_export_headers(tmp_path: Path):
    source = tmp_path / "fish.csv"
    pd.DataFrame(
        {"accession": ["ISO-01"], "arg": ["blaTEM-1"], "pid": ["99.1"], "cov": ["95"]}
    ).to_csv(source, index=False)
    mapping = tmp_path / "mapping.json"
    mapping.write_text(
        json.dumps(
            {
                "columns": {
                    "sample_id": "accession",
                    "gene": "arg",
                    "identity_pct": "pid",
                    "coverage_pct": "cov",
                },
                "defaults": {"source_system": "FishAMR-Link"},
            }
        )
    )
    target = tmp_path / "normalized.csv"
    result = normalize_fishamr_calls(source, target, mapping)
    assert result.loc[0, "sample_id"] == "ISO-01"
    assert result.loc[0, "gene"] == "blaTEM-1"
    assert result.loc[0, "source_system"] == "FishAMR-Link"
    assert target.exists()


def test_epiresnet_adapter_keeps_prediction_type_separate(tmp_path: Path):
    source = tmp_path / "epi.csv"
    pd.DataFrame({"id": ["ISO-01"], "drug": ["Cefotaxime"], "mic_hat": ["4"]}).to_csv(
        source, index=False
    )
    mapping = tmp_path / "mapping.json"
    mapping.write_text(
        json.dumps(
            {
                "columns": {"sample_id": "id", "antimicrobial": "drug", "predicted_mic": "mic_hat"},
                "defaults": {"model_name": "EpiResNet"},
            }
        )
    )
    target = tmp_path / "predictions.csv"
    result = normalize_epiresnet_predictions(source, target, mapping)
    assert result.loc[0, "predicted_mic"] == "4"
    assert "result" not in result.columns
    assert result.loc[0, "model_name"] == "EpiResNet"


def test_adapter_fails_on_missing_mapped_source_column(tmp_path: Path):
    source = tmp_path / "fish.csv"
    pd.DataFrame({"id": ["ISO-01"]}).to_csv(source, index=False)
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"columns": {"sample_id": "id", "gene": "missing_gene"}}))
    with pytest.raises(ConfigurationError):
        normalize_fishamr_calls(source, tmp_path / "out.csv", mapping)


def test_csv_reader_rejects_duplicate_headers_before_pandas_mangles_them(tmp_path: Path):
    from amr_surveil.exceptions import InputFileError
    from amr_surveil.io import read_csv

    source = tmp_path / "duplicate_header.csv"
    source.write_text("sample_id, sample_id \nISO-1,ISO-2\n", encoding="utf-8")
    with pytest.raises(InputFileError, match="duplicate column"):
        read_csv(source, table="duplicate_test")


def test_csv_reader_rejects_blank_headers_and_empty_files(tmp_path: Path):
    from amr_surveil.exceptions import InputFileError
    from amr_surveil.io import read_csv

    blank = tmp_path / "blank_header.csv"
    blank.write_text("sample_id,,organism\nA,,E. coli\n", encoding="utf-8")
    with pytest.raises(InputFileError, match="blank column"):
        read_csv(blank)

    empty = tmp_path / "empty.csv"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(InputFileError, match="no header"):
        read_csv(empty)
