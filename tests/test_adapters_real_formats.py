"""Adapters against fixtures that follow the documented export layouts of the upstream tools.

The fixtures are hand-written from the published AMRFinderPlus column headers, the FishAMR-Link
source (table names, columns) and the EpiResNet training-script outputs. They are NOT captured
from live tool runs; see docs/LIMITATIONS.md.
"""

import json
import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from amr_surveil.adapters.amrfinder import normalize_amrfinder_reports
from amr_surveil.adapters.epiresnet import normalize_epiresnet_test_predictions
from amr_surveil.adapters.fishamr import (
    infer_fishamr_tier,
    normalize_fishamr_feature_table,
    read_fishamr_warehouse,
)
from amr_surveil.exceptions import ConfigurationError, InputFileError

V3 = [
    "Name",
    "Protein identifier",
    "Contig id",
    "Start",
    "Stop",
    "Strand",
    "Gene symbol",
    "Sequence name",
    "Scope",
    "Element type",
    "Element subtype",
    "Class",
    "Subclass",
    "Method",
    "Target length",
    "Reference sequence length",
    "% Coverage of reference sequence",
    "% Identity to reference sequence",
    "Alignment length",
    "Accession of closest sequence",
    "Name of closest sequence",
    "HMM id",
    "HMM description",
]


def tsv(path: Path, rows: list[dict], header=V3):
    lines = ["\t".join(header)] + ["\t".join(str(r.get(h, "NA")) for h in header) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def v3_row(
    name,
    gene,
    etype="AMR",
    method="EXACTX",
    ident="100.00",
    cov="100.00",
    cls="BETA-LACTAM",
    sub="CEPHALOSPORIN",
):
    return {
        "Name": name,
        "Contig id": "c1",
        "Start": 10,
        "Stop": 900,
        "Gene symbol": gene,
        "Sequence name": gene + " enzyme",
        "Element type": etype,
        "Class": cls,
        "Subclass": sub,
        "Method": method,
        "% Coverage of reference sequence": cov,
        "% Identity to reference sequence": ident,
    }


def test_amrfinder_v3_basic_and_element_type_filter(tmp_path):
    f = tsv(
        tmp_path / "s1.tsv",
        [
            v3_row("S1", "blaCTX-M-15"),
            v3_row("S1", "acrF", etype="STRESS"),
            v3_row("S1", "tetA", method="PARTIALX", ident="97.1"),
        ],
    )
    summary = normalize_amrfinder_reports(
        [f],
        tmp_path / "calls.csv",
        tmp_path / "scr.csv",
        tool_version="3.12.8",
        database_version="2024-01-31.1",
    )
    calls = pd.read_csv(tmp_path / "calls.csv")
    assert set(calls["gene"]) == {"blaCTX-M-15", "tetA"}
    assert set(calls["evidence_tier"]) == {"curated_database"}
    assert calls.loc[calls.gene == "tetA", "hit_method"].iloc[0] == "PARTIALX"
    assert summary["excluded_by_element_type"] == {"STRESS": 1}
    scr = pd.read_csv(tmp_path / "scr.csv")
    assert (
        scr.loc[0, "tool_version"] == "3.12.8" and scr.loc[0, "database_version"] == "2024-01-31.1"
    )


def test_amrfinder_v4_headers_accepted(tmp_path):
    header = [
        "Name",
        "Contig id",
        "Start",
        "Stop",
        "Element symbol",
        "Element name",
        "Type",
        "Subtype",
        "Class",
        "Subclass",
        "Method",
        "% Coverage of reference",
        "% Identity to reference",
    ]
    f = tmp_path / "x.tsv"
    f.write_text(
        "\t".join(header)
        + "\n"
        + "\t".join(
            [
                "S9",
                "c",
                "1",
                "2",
                "qnrS1",
                "quinolone resistance",
                "AMR",
                "AMR",
                "QUINOLONE",
                "QUINOLONE",
                "EXACTX",
                "100",
                "100",
            ]
        )
        + "\n"
    )
    normalize_amrfinder_reports([f], tmp_path / "c.csv")
    calls = pd.read_csv(tmp_path / "c.csv")
    assert calls.loc[0, "gene"] == "qnrS1" and calls.loc[0, "sample_id"] == "S9"


def test_header_only_report_is_a_screened_negative_not_missing(tmp_path):
    f = tsv(tmp_path / "clean_genome.tsv", [])
    normalize_amrfinder_reports([f], tmp_path / "c.csv", tmp_path / "s.csv")
    assert len(pd.read_csv(tmp_path / "c.csv")) == 0
    assert pd.read_csv(tmp_path / "s.csv")["sample_id"].tolist() == ["clean_genome"]


def test_zero_byte_report_is_an_error_not_a_negative(tmp_path):
    f = tmp_path / "broken.tsv"
    f.write_text("")
    with pytest.raises(InputFileError, match="empty"):
        normalize_amrfinder_reports([f], tmp_path / "c.csv")


def test_non_amrfinder_file_and_duplicate_sample_rejected(tmp_path):
    f = tmp_path / "random.tsv"
    f.write_text("a\tb\n1\t2\n")
    with pytest.raises(InputFileError, match="gene column"):
        normalize_amrfinder_reports([f], tmp_path / "c.csv")
    a = tsv(tmp_path / "a.tsv", [v3_row("DUP", "blaTEM-1")])
    b = tsv(tmp_path / "b.tsv", [v3_row("DUP", "tetA")])
    with pytest.raises(ConfigurationError, match="more than one report"):
        normalize_amrfinder_reports([a, b], tmp_path / "c.csv")


def test_sample_ids_argument_used_when_no_name_column(tmp_path):
    header = [h for h in V3 if h != "Name"]
    f = tsv(tmp_path / "r.tsv", [v3_row("", "blaTEM-1")], header=header)
    normalize_amrfinder_reports([f], tmp_path / "c.csv", sample_ids=["ISO-7"])
    assert pd.read_csv(tmp_path / "c.csv").loc[0, "sample_id"] == "ISO-7"
    with pytest.raises(ConfigurationError):
        normalize_amrfinder_reports([f], tmp_path / "c.csv", sample_ids=["a", "b"])


# ------------------------------------------------------------------ FishAMR-Link
@pytest.mark.parametrize(
    "klass,mech,tier",
    [
        ("BETA-LACTAM", "EXACTX", "curated_database"),
        ("QUINOLONE", "POINTX", "curated_database"),
        ("Homology Match", "Identified via remote blastn (98% id)", "homology_search"),
        ("Aminoglycoside", "Enzymatic Inactivation", "heuristic_motif"),
        ("X", "something new", "unspecified"),
        (None, None, "unspecified"),
    ],
)
def test_fishamr_tier_inference(klass, mech, tier):
    assert infer_fishamr_tier("g", klass, mech) == tier


def test_fishamr_feature_table(tmp_path):
    src = tmp_path / "amr_genes.csv"
    pd.DataFrame(
        [
            {
                "gene": "blaCTX-M-15",
                "class": "BETA-LACTAM",
                "mechanism": "EXACTX",
                "feature_type": "AMR_Gene",
                "confidence": 1.0,
            },
            {
                "gene": "sul1",
                "class": "Sulfonamide",
                "mechanism": "Dihydropteroate Synthase Bypass",
                "feature_type": "AMR_Gene",
                "confidence": 0.89,
            },
            {
                "gene": "IncFIB",
                "class": "Plasmid",
                "mechanism": "Autonomous Replication Initiation",
                "feature_type": "MGE",
                "confidence": 0.92,
            },
            {
                "gene": "IS26",
                "class": "Transposon",
                "mechanism": "Replicative Transposition",
                "feature_type": "MGE",
                "confidence": 0.92,
            },
            {
                "gene": "iucA",
                "class": "Virulence",
                "mechanism": "Iron Acquisition Sluice",
                "feature_type": "Virulence",
                "confidence": 0.9,
            },
        ]
    ).to_csv(src, index=False)
    with pytest.raises(ConfigurationError, match="sample_id"):
        normalize_fishamr_feature_table(src, tmp_path / "c.csv")
    summary = normalize_fishamr_feature_table(
        src, tmp_path / "c.csv", tmp_path / "p.csv", sample_id="ISO-1", organism="Vibrio"
    )
    calls = pd.read_csv(tmp_path / "c.csv")
    assert dict(zip(calls.gene, calls.evidence_tier, strict=True)) == {
        "blaCTX-M-15": "curated_database",
        "sul1": "heuristic_motif",
    }
    assert "confidence" not in calls.columns  # fixed constants are not statistical confidence
    assert pd.read_csv(tmp_path / "p.csv")["replicon"].tolist() == [
        "IncFIB"
    ]  # IS26 is not a replicon
    assert summary["excluded_feature_types"] == {"Virulence": 1, "MGE": 1}


def make_warehouse(path: Path):
    con = sqlite3.connect(path)
    con.executescript(
        """CREATE TABLE isolates(id INTEGER PRIMARY KEY, accession TEXT, organism TEXT, host_species TEXT, country TEXT,
              collection_year INTEGER, isolation_source TEXT, latitude REAL, longitude REAL);
           CREATE TABLE genomic_sequences(id INTEGER PRIMARY KEY, isolate_id INTEGER, detected_amr_genes TEXT);"""
    )
    con.executemany(
        "INSERT INTO isolates VALUES (?,?,?,?,?,?,?,?,?)",
        [
            (
                1,
                "NZ_CP012345.1",
                "Escherichia coli",
                "fish",
                "X",
                2020,
                "gut",
                1.0,
                2.0,
            ),  # seeded simulated baseline
            (2, "REAL-001", "Vibrio parahaemolyticus", "shrimp", "Y", 2023, "pond", 3.0, 4.0),
        ],
    )
    con.executemany(
        "INSERT INTO genomic_sequences VALUES (?,?,?)",
        [
            (1, 1, json.dumps([{"gene": "tetA"}])),
            (2, 2, json.dumps([{"gene": "blaCARB"}, "qnrVC"])),
        ],
    )
    con.commit()
    con.close()


def test_warehouse_export_skips_seeded_baseline_and_coordinates(tmp_path):
    db = tmp_path / "w.db"
    make_warehouse(db)
    s = read_fishamr_warehouse(db, tmp_path / "s.csv", tmp_path / "c.csv")
    samples, calls = pd.read_csv(tmp_path / "s.csv"), pd.read_csv(tmp_path / "c.csv")
    assert (
        samples["sample_id"].tolist() == ["REAL-001"]
        and s["n_seeded_baseline_isolates_skipped"] == 1
    )
    assert sorted(calls["gene"]) == ["blaCARB", "qnrVC"]
    assert not {"latitude", "longitude"} & set(samples.columns)
    s2 = read_fishamr_warehouse(
        db, tmp_path / "s2.csv", tmp_path / "c2.csv", include_seeded_baseline=True
    )
    assert s2["n_isolates_exported"] == 2
    assert (
        pd.read_csv(tmp_path / "s2.csv")
        .query("sample_id == 'NZ_CP012345.1'")["data_license"]
        .iloc[0]
        == "synthetic"
    )


def test_warehouse_is_opened_read_only_and_validated(tmp_path):
    db = tmp_path / "w.db"
    make_warehouse(db)
    before = db.read_bytes()
    read_fishamr_warehouse(db, tmp_path / "s.csv", tmp_path / "c.csv")
    assert db.read_bytes() == before
    other = tmp_path / "other.db"
    con = sqlite3.connect(other)
    con.executescript("CREATE TABLE t(a)")
    con.close()
    with pytest.raises(InputFileError, match="expected tables"):
        read_fishamr_warehouse(other, tmp_path / "s.csv", tmp_path / "c.csv")
    with pytest.raises(InputFileError, match="not found"):
        read_fishamr_warehouse(tmp_path / "nope.db", tmp_path / "s.csv", tmp_path / "c.csv")


# ------------------------------------------------------------------- EpiResNet
def test_epiresnet_test_predictions(tmp_path):
    csv = tmp_path / "test_predictions.csv"
    pd.DataFrame(
        {
            "isolate_id": ["A", "A", "B", "C"],
            "antibiotic_name": ["cefotaxime"] * 4,
            "protein_cluster": ["p1", "p2", "p1", "p3"],
            "binary_valid": [True, True, True, False],
            "probability_resistant": [0.9, 0.2, 0.7, None],
            "prediction_resistant": [1, 0, 1, None],
        }
    ).to_csv(csv, index=False)
    (tmp_path / "final_metrics.json").write_text(
        json.dumps(
            {
                "selection": {
                    "selected_threshold_on_validation": 0.41,
                    "probability_temperature": 1.3,
                }
            }
        )
    )
    (tmp_path / "run_manifest.json").write_text(
        json.dumps({"input_sha256": "abc123", "config": {"split_strategy": "cluster"}})
    )
    ckpt = tmp_path / "model.pt"
    ckpt.write_bytes(b"weights")
    summary = normalize_epiresnet_test_predictions(
        csv,
        tmp_path / "out.csv",
        final_metrics_path=tmp_path / "final_metrics.json",
        run_manifest_path=tmp_path / "run_manifest.json",
        checkpoint_path=ckpt,
    )
    out = pd.read_csv(tmp_path / "out.csv")
    assert len(out) == 3 and summary["n_rows_without_model_output"] == 1
    assert set(out["predicted_class"]) == {"predicted_resistant", "predicted_susceptible"}
    assert out["evaluation_split"].unique().tolist() == ["test"]
    assert out["decision_threshold"].unique().tolist() == [0.41] and out[
        "probability_temperature"
    ].unique().tolist() == [1.3]
    assert out["split_strategy"].unique().tolist() == ["cluster"] and out[
        "training_data_version"
    ].unique().tolist() == ["abc123"]
    assert out["checkpoint_id"].nunique() == 1 and len(out["checkpoint_id"].iloc[0]) == 64
    assert (out["prediction_unit"] == "protein-antibiotic pair").all()
    assert out["predicted_mic"].isna().all()  # EpiResNet v5 exports no MIC; none is invented


def test_epiresnet_missing_columns_rejected(tmp_path):
    csv = tmp_path / "x.csv"
    pd.DataFrame({"isolate_id": ["A"]}).to_csv(csv, index=False)
    with pytest.raises(InputFileError, match="missing EpiResNet"):
        normalize_epiresnet_test_predictions(csv, tmp_path / "o.csv")
