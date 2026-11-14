import io
import sys
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pandas as pd
import pytest

from amr_surveil.exceptions import ConfigurationError
from amr_surveil.synthetic import generate_synthetic_dataset
from amr_surveil.ui_support import (
    UPLOAD_SLOTS,
    build_config,
    error_bars,
    read_result,
    run_demo,
    run_from_files,
    zip_outputs,
)

ROOT = Path(__file__).resolve().parents[1]


def test_build_config_defaults_and_validation():
    cfg = build_config()
    assert cfg["quality"]["strict"] and set(cfg["analysis"]["excluded_evidence_tiers"]) == {
        "heuristic_motif",
        "homology_search",
    }
    assert build_config(excluded_evidence_tiers=[])["analysis"]["excluded_evidence_tiers"] == []
    with pytest.raises(ConfigurationError):
        build_config(deduplication_rule="nope")
    with pytest.raises(ConfigurationError):
        build_config(excluded_evidence_tiers=["bogus"])


def test_upload_slots_match_pipeline_inputs():
    assert [k for k, _, req in UPLOAD_SLOTS if req] == ["samples"]
    assert len({k for k, _, _ in UPLOAD_SLOTS}) == len(UPLOAD_SLOTS) == 11


def test_run_demo_and_zip(tmp_path):
    summary, out = run_demo(tmp_path, n_samples=40)
    assert summary["n_samples"] == 40 and (out / "report.md").is_file()
    names = zipfile.ZipFile(io.BytesIO(zip_outputs(out))).namelist()
    assert "report.md" in names and "run_manifest.json" in names


def test_run_from_files_uses_ui_config(tmp_path):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=40)
    cfg = build_config(sample_selection_documented=True, selection_description="all")
    summary = run_from_files(paths, cfg, tmp_path / "w" / "results")
    assert summary["genotype_phenotype_decision"] == "GO"


def test_read_result_tolerates_missing_and_headerless(tmp_path):
    assert read_result(tmp_path, "nope.csv").empty
    (tmp_path / "blank.csv").write_text("")
    assert read_result(tmp_path, "blank.csv").empty


def test_error_bars_never_negative():
    f = pd.DataFrame(
        {
            "resistant_pct": [50.0, 0.0],
            "wilson_low_pct": [40.0, 0.0],
            "wilson_high_pct": [60.0, 30.0],
        }
    )
    up, down = error_bars(f)
    assert (up >= 0).all() and (down >= 0).all() and up.tolist() == [10.0, 30.0]


# ---------------------------------------------------------------- Streamlit script smoke test
def fake_streamlit(uploads=None, click=None):
    st = MagicMock(name="streamlit")
    st.session_state = {}
    # Column objects delegate to the same stub so ``col.selectbox(...)`` returns real values.
    st.columns.side_effect = lambda spec: [
        st for _ in range(spec if isinstance(spec, int) else len(spec))
    ]
    st.number_input.side_effect = lambda label, min_value=None, max_value=None, value=None, **kw: (
        value
    )
    st.radio.side_effect = lambda label, options, **kw: options[0]
    st.tabs.side_effect = lambda labels: [MagicMock() for _ in labels]
    st.file_uploader.side_effect = lambda label, **kw: (uploads or {}).get(
        kw["key"].removeprefix("up_")
    )
    st.button.side_effect = lambda label, **kw: label == click
    st.checkbox.side_effect = lambda label, value=False, **kw: value
    st.selectbox.side_effect = lambda label, options, index=0, **kw: options[index]
    st.multiselect.side_effect = lambda label, options, default=None, **kw: list(default or [])
    st.text_input.return_value = ""
    st.text_area.return_value = ""
    st.stop.side_effect = SystemExit
    return st


@pytest.fixture
def stubbed(monkeypatch):
    px = MagicMock(name="plotly.express")
    monkeypatch.setitem(sys.modules, "plotly", SimpleNamespace(express=px))
    monkeypatch.setitem(sys.modules, "plotly.express", px)

    def install(st):
        errors = SimpleNamespace(
            StreamlitAPIException=type("StreamlitAPIException", (Exception,), {})
        )
        st.errors = errors
        monkeypatch.setitem(sys.modules, "streamlit", st)
        monkeypatch.setitem(sys.modules, "streamlit.errors", errors)

    return install


def run_app():
    import runpy

    runpy.run_path(str(ROOT / "app.py"), run_name="__main__")


def test_app_script_renders_demo_without_errors(stubbed):
    st = fake_streamlit()
    stubbed(st)
    run_app()
    st.error.assert_not_called()
    assert st.dataframe.call_count >= 10 and st.download_button.called
    assert st.session_state["amr_surveil_summary"]["n_samples"] == 120
    px = sys.modules["plotly.express"]
    assert px.scatter.called  # AST indicators chart
    assert px.bar.called  # AST re-check chart (Plotly, not altair-backed st.bar_chart)
    assert any(
        c.args and c.args[0] == "Download antibiogram CSV"
        for c in st.download_button.call_args_list
    )


def test_app_script_handles_upload_run_and_error_path(stubbed, tmp_path):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=40)
    up = {k: SimpleNamespace(getvalue=p.read_bytes) for k, p in paths.items()}
    st = fake_streamlit(uploads=up, click="Run analysis")
    stubbed(st)
    run_app()
    st.error.assert_not_called()
    assert st.session_state["amr_surveil_summary"]["n_samples"] == 40

    bad = generate_synthetic_dataset(tmp_path / "bad", n_samples=40, inject_defects=True)
    up = {
        "samples": SimpleNamespace(getvalue=bad["samples"].read_bytes),
        "ast": SimpleNamespace(getvalue=bad["ast"].read_bytes),
    }
    st = fake_streamlit(uploads=up, click="Run analysis")
    stubbed(st)
    with pytest.raises(SystemExit):  # st.stop(): nothing to render after a failed strict run
        run_app()
    st.error.assert_called()  # strict gate message surfaced, no crash
    assert "Strict quality gate" in st.error.call_args.args[0]


def test_read_table_upload_sniffs_delimiter_and_keeps_text():
    from amr_surveil.ui_support import read_table_upload

    tsv = read_table_upload(b"a\tb\n001\tNA\n", "x.tsv")
    assert list(tsv.columns) == ["a", "b"] and tsv.loc[0, "a"] == "001" and tsv.loc[0, "b"] == "NA"
    csv = read_table_upload("﻿a,b\n1,2\n".encode(), "x.csv")
    assert list(csv.columns) == ["a", "b"]
    assert read_table_upload(b"", "e.csv").empty


def test_context_helpers():
    from amr_surveil.ui_support import context_labels, filter_context

    frame = pd.DataFrame(
        {
            "method": ["mic", "disk"],
            "breakpoint_standard": ["CLSI", "EUCAST"],
            "breakpoint_version": ["2019", "2020"],
        }
    )
    labels = context_labels(frame)
    assert labels == ["disk | EUCAST | 2020", "mic | CLSI | 2019"]
    assert (
        len(filter_context(frame, "mic | CLSI | 2019")) == 1
        and len(filter_context(frame, None)) == 2
    )
    assert context_labels(pd.DataFrame()) == []


def test_analyse_public_dataset_end_to_end(tmp_path):
    from amr_surveil.datasources import parse_bvbrc
    from amr_surveil.ui_support import analyse_public_dataset, build_config, public_inputs_zip

    amr = [
        {
            "genome_id": f"1.{i}",
            "genome_name": "Escherichia coli x",
            "antibiotic": "ampicillin",
            "resistant_phenotype": "Resistant" if i % 2 else "Susceptible",
            "laboratory_typing_method": "MIC",
            "testing_standard": "CLSI",
            "testing_standard_year": "2019",
        }
        for i in range(10)
    ]
    genome = [
        {
            "genome_id": f"1.{i}",
            "species": "Escherichia coli",
            "isolation_country": "Kenya",
            "collection_year": "2019",
            "host_name": "Homo sapiens",
        }
        for i in range(10)
    ]
    ds = parse_bvbrc(amr, genome, organism="Escherichia coli")
    summary, out = analyse_public_dataset(ds, tmp_path, build_config(strict=False))
    assert summary["n_ast_records"] == 10 and (out / "resistance_by_group.csv").is_file()
    assert (
        "public_data_provenance.json"
        in zipfile.ZipFile(io.BytesIO(public_inputs_zip(ds, tmp_path / "z"))).namelist()
    )
    import pytest as _pt

    from amr_surveil.datasources import parse_ecdc_atlas

    agg_only = parse_ecdc_atlas(pd.DataFrame({"Indicator": [], "Time": [], "NumValue": []}))
    with _pt.raises(ValueError):
        analyse_public_dataset(agg_only, tmp_path / "agg", build_config(strict=False))


def _public_dataset():
    from amr_surveil.datasources import parse_bvbrc

    amr = [
        {
            "genome_id": f"1.{i}",
            "genome_name": "Escherichia coli x",
            "antibiotic": "ampicillin",
            "resistant_phenotype": "Resistant" if i % 3 else "Susceptible",
            "laboratory_typing_method": "MIC",
            "testing_standard": "CLSI",
            "testing_standard_year": "2019",
        }
        for i in range(30)
    ]
    genome = [
        {
            "genome_id": f"1.{i}",
            "species": "Escherichia coli",
            "isolation_country": "Kenya",
            "collection_year": "2019",
            "host_name": "Homo sapiens",
        }
        for i in range(30)
    ]
    return parse_bvbrc(amr, genome, organism="Escherichia coli")


def test_app_hub_retrieve_and_analyse_flow(stubbed, monkeypatch):
    import amr_surveil.datasources as ds_pkg

    monkeypatch.setattr(ds_pkg, "fetch_bvbrc", lambda *a, **k: _public_dataset())
    st = fake_streamlit(click="Retrieve from BV-BRC")
    stubbed(st)
    run_app()
    st.error.assert_not_called()
    held = st.session_state["hub_dataset"]
    assert len(held.samples) == 30

    # second script run: the held dataset is analysed and the results replace the demo
    st2 = fake_streamlit(click="Analyse this dataset")
    st2.session_state["hub_dataset"] = held
    stubbed(st2)
    run_app()
    st2.error.assert_not_called()
    assert st2.rerun.called and st2.session_state["amr_surveil_summary"]["n_ast_records"] == 30


def test_app_hub_reports_fetch_errors_without_crashing(stubbed, monkeypatch):
    import amr_surveil.datasources as ds_pkg
    from amr_surveil.datasources import FetchError

    def boom(*a, **k):
        raise FetchError("could not reach www.bv-brc.org")

    monkeypatch.setattr(ds_pkg, "fetch_bvbrc", boom)
    st = fake_streamlit(click="Retrieve from BV-BRC")
    stubbed(st)
    run_app()
    st.error.assert_called_once()
    assert "could not reach" in st.error.call_args.args[0]


def test_app_survives_undated_public_dataset(stubbed, tmp_path):
    from amr_surveil.datasources import parse_bvbrc_export
    from amr_surveil.ui_support import analyse_public_dataset

    table = pd.DataFrame(
        {
            "Genome ID": [f"1.{i}" for i in range(12)],
            "Genome Name": "Escherichia coli x",
            "Antibiotic": "ampicillin",
            "Resistant Phenotype": ["Resistant", "Susceptible"] * 6,
        }
    )
    _, out = analyse_public_dataset(parse_bvbrc_export(table), tmp_path, build_config(strict=False))
    st = fake_streamlit()
    st.session_state["amr_surveil_result_dir"] = str(out)
    st.session_state["amr_surveil_summary"] = {"n_samples": 12}
    stubbed(st)
    run_app()
    st.error.assert_not_called()
