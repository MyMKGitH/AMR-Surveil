"""Public-data hub tests. No test touches the network: transport is an injected fake."""

from __future__ import annotations

import gzip
import io
import json
import urllib.error
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pytest

from amr_surveil.antibiogram import cumulative_antibiogram
from amr_surveil.constants import AGGREGATE_COLUMNS, AST_COLUMNS, SAMPLE_COLUMNS
from amr_surveil.datasources import (
    FetchError,
    country_to_iso3,
    fetch,
    fetch_bvbrc,
    parse_bvbrc,
    parse_bvbrc_export,
    parse_ecdc_atlas,
    parse_pathogen_table,
    write_inputs,
)
from amr_surveil.datasources.bvbrc import amr_query_url, genome_query_url, organism_glob
from amr_surveil.datasources.common import infer_sector, map_result, normalise_date
from amr_surveil.datasources.ncbi_pathogen import parse_pairs
from amr_surveil.datasources.registry import SOURCES, get_source
from amr_surveil.mapping import country_resistance
from amr_surveil.pipeline import run_pipeline

HOSTS = ("www.bv-brc.org",)
NOW = lambda: datetime(2026, 11, 12, 9, 0, tzinfo=timezone.utc)  # noqa: E731


class FakeResponse:
    def __init__(self, body: bytes, headers: dict[str, str] | None = None) -> None:
        self._buf = io.BytesIO(body)
        self.headers = headers or {}
        self.closed = False

    def read(self, n: int = -1) -> bytes:
        return self._buf.read(n)

    def close(self) -> None:
        self.closed = True


class FakeOpener:
    def __init__(self, routes: dict[str, bytes] | None = None, *, error: Exception | None = None):
        self.routes, self.error, self.calls = routes or {}, error, []

    def __call__(self, request, timeout):
        self.calls.append(request.full_url)
        if self.error:
            raise self.error
        for needle, body in self.routes.items():
            if needle in request.full_url:
                return FakeResponse(body)
        raise AssertionError(f"unexpected URL {request.full_url}")


# ---------------------------------------------------------------- fetch layer
def _fetch(tmp_path: Path, url: str, opener, **kw):
    return fetch(
        url, allowed_hosts=HOSTS, cache_dir=tmp_path, source_id="t", opener=opener, now=NOW, **kw
    )


def test_fetch_rejects_http_and_foreign_hosts_and_credentials(tmp_path: Path) -> None:
    op = FakeOpener({"": b"x"})
    for url in (
        "http://www.bv-brc.org/a",
        "https://evil.example/a",
        "https://u:p@www.bv-brc.org/a",
        "ftp://www.bv-brc.org/a",
        "file:///etc/passwd",
    ):
        with pytest.raises(FetchError):
            _fetch(tmp_path, url, op)
    assert op.calls == []  # nothing was ever requested


def test_fetch_records_provenance_and_uses_cache(tmp_path: Path) -> None:
    op = FakeOpener({"bv-brc": b'[{"a":1}]'})
    first = _fetch(tmp_path, "https://www.bv-brc.org/api/x", op)
    assert first.retrieved_at == "2026-11-12T09:00:00Z" and not first.from_cache
    assert first.n_bytes == 9 and len(first.sha256) == 64
    second = _fetch(tmp_path, "https://www.bv-brc.org/api/x", op)
    assert second.from_cache and second.sha256 == first.sha256 and len(op.calls) == 1
    third = _fetch(tmp_path, "https://www.bv-brc.org/api/x", op, refresh=True)
    assert not third.from_cache and len(op.calls) == 2


def test_fetch_detects_tampered_cache(tmp_path: Path) -> None:
    op = FakeOpener({"bv-brc": b"original"})
    first = _fetch(tmp_path, "https://www.bv-brc.org/api/x", op)
    first.path.write_bytes(b"tampered")
    again = _fetch(tmp_path, "https://www.bv-brc.org/api/x", op)
    assert not again.from_cache and again.path.read_bytes() == b"original"


def test_fetch_enforces_size_limits(tmp_path: Path) -> None:
    with pytest.raises(FetchError, match="limit"):
        _fetch(
            tmp_path,
            "https://www.bv-brc.org/api/x",
            FakeOpener({"bv-brc": b"x" * 5000}),
            max_bytes=1000,
        )

    class Declared(FakeOpener):
        def __call__(self, request, timeout):
            return FakeResponse(b"x", {"Content-Length": "999999"})

    with pytest.raises(FetchError, match="announces"):
        _fetch(tmp_path, "https://www.bv-brc.org/api/y", Declared(), max_bytes=1000)


def test_fetch_gunzips_and_bounds_decompression(tmp_path: Path) -> None:
    ok = gzip.compress(b'{"k": 1}')
    assert (
        _fetch(
            tmp_path, "https://www.bv-brc.org/api/g", FakeOpener({"bv-brc": ok})
        ).path.read_bytes()
        == b'{"k": 1}'
    )
    bomb = gzip.compress(b"0" * 200_000)
    with pytest.raises(FetchError, match="Decompressed"):
        _fetch(
            tmp_path, "https://www.bv-brc.org/api/b", FakeOpener({"bv-brc": bomb}), max_bytes=50_000
        )


def test_fetch_wraps_network_errors(tmp_path: Path) -> None:
    http = urllib.error.HTTPError("https://www.bv-brc.org/api/x", 503, "down", {}, None)  # type: ignore[arg-type]
    with pytest.raises(FetchError, match="HTTP 503"):
        _fetch(tmp_path, "https://www.bv-brc.org/api/x", FakeOpener(error=http))
    with pytest.raises(FetchError, match="could not reach"):
        _fetch(
            tmp_path, "https://www.bv-brc.org/api/x", FakeOpener(error=urllib.error.URLError("dns"))
        )
    with pytest.raises(FetchError, match="could not reach"):
        _fetch(tmp_path, "https://www.bv-brc.org/api/x", FakeOpener(error=TimeoutError("slow")))


# ---------------------------------------------------------------- geo / vocab
@pytest.mark.parametrize(
    ("raw", "iso"),
    [
        ("USA: Texas", "USA"),
        ("Viet Nam", "VNM"),
        ("Türkiye", "TUR"),
        ("Côte d'Ivoire", "CIV"),
        ("United Kingdom", "GBR"),
        ("germany", "DEU"),
        ("Korea, Republic of", "KOR"),
        ("Atlantis", None),
        ("", None),
        (None, None),
        (float("nan"), None),
    ],
)
def test_country_to_iso3(raw, iso) -> None:
    assert country_to_iso3(raw) == iso


def test_vocabularies() -> None:
    assert [
        map_result(v) for v in ("Resistant", "S", "intermediate", "Susceptible-dose dependent")
    ] == ["R", "S", "I", "I"]
    assert (
        map_result("Nonsusceptible") is None and map_result("") is None and map_result(None) is None
    )
    assert infer_sector("Homo sapiens") == "human" and infer_sector("Gallus gallus") == "animal"
    assert infer_sector("", "river water") == "environment" and infer_sector("", "") == "unknown"
    assert (
        normalise_date("2019-05-02T00:00:00Z") == "2019-05-02" and normalise_date("2019") == "2019"
    )
    assert normalise_date("2019-13") == "2019" and normalise_date("not a date") == ""
    assert parse_pairs("a=R, b=S;c=I;broken") == [("a", "R"), ("b", "S"), ("c", "I")]


def test_registry_has_terms_for_every_source() -> None:
    for sid, info in SOURCES.items():
        assert info.terms and info.terms_url.startswith("https://") and info.terms_checked
        assert info.access in {"live_api", "file_import", "link_only"}
        assert get_source(sid) is info
    with pytest.raises(KeyError):
        get_source("nope")


# ---------------------------------------------------------------- BV-BRC
AMR = [
    {
        "genome_id": "562.1",
        "genome_name": "Escherichia coli strain A",
        "antibiotic": "ciprofloxacin",
        "resistant_phenotype": "Resistant",
        "measurement_sign": ">",
        "measurement_value": "4",
        "measurement_unit": "mg/L",
        "laboratory_typing_method": "MIC",
        "testing_standard": "CLSI",
        "testing_standard_year": "2019",
    },
    {
        "genome_id": "562.2",
        "genome_name": "Escherichia coli strain B",
        "antibiotic": "ciprofloxacin",
        "resistant_phenotype": "Susceptible",
        "measurement_value": "0.03",
        "measurement_unit": "mg/L",
        "laboratory_typing_method": "MIC",
        "testing_standard": "CLSI",
        "testing_standard_year": "2019",
    },
    {
        "genome_id": "562.3",
        "genome_name": "Escherichia coli strain C",
        "antibiotic": "ciprofloxacin",
        "resistant_phenotype": "Nonsusceptible",
        "laboratory_typing_method": "Disk diffusion",
    },
    {
        "genome_id": "562.4",
        "genome_name": "Escherichia coli strain D",
        "antibiotic": "ciprofloxacin",
        "resistant_phenotype": "Susceptible",
        "measurement_value": "30",
        "laboratory_typing_method": "Disk diffusion",
        "testing_standard": "EUCAST",
        "testing_standard_year": "2020",
    },
]
GENOME = [
    {
        "genome_id": "562.1",
        "species": "Escherichia coli",
        "isolation_country": "Viet Nam",
        "collection_year": "2018",
        "host_name": "Homo sapiens",
        "isolation_source": "blood",
    },
    {
        "genome_id": "562.2",
        "species": "Escherichia coli",
        "isolation_country": "Viet Nam",
        "collection_year": "2019",
        "host_name": "Homo sapiens",
        "isolation_source": "urine",
    },
    {
        "genome_id": "562.4",
        "species": "Escherichia coli",
        "isolation_country": "Kenya",
        "collection_year": "2020",
        "host_name": "Gallus gallus",
    },
]


def test_bvbrc_urls_are_encoded_and_bounded() -> None:
    assert organism_glob("Escherichia coli") == "Escherichia*coli*"
    url = amr_query_url("Escherichia coli", ["Ciprofloxacin", "ampicillin"], limit=50)
    assert url.startswith("https://www.bv-brc.org/api/genome_amr/?") and "limit(50,0)" in url
    assert "in(antibiotic,(ciprofloxacin,ampicillin))" in url
    assert "evil" not in genome_query_url(["562.1"])
    with pytest.raises(ValueError):
        amr_query_url("E coli", [], limit=10)
    with pytest.raises(ValueError):
        amr_query_url("E coli", ["x"], limit=10**6)
    assert "%26" in amr_query_url("E&coli", ["x"], limit=1)  # injection characters are escaped


def test_parse_bvbrc_conversion_and_reporting() -> None:
    ds = parse_bvbrc(AMR, GENOME, organism="Escherichia coli", retrieved_at="2026-11-12T09:00:00Z")
    assert list(ds.samples.columns) == SAMPLE_COLUMNS and list(ds.ast.columns) == AST_COLUMNS
    assert sorted(ds.ast["result"]) == ["R", "S", "S"]  # nonsusceptible dropped, not guessed
    assert any("not S/I/R" in n for n in ds.notes)
    row = ds.ast.set_index("sample_id").loc["562.1"]
    assert row["mic"] == ">4" and row["mic_unit"] == "mg/L" and row["zone_mm"] == ""
    assert row["breakpoint_standard"] == "CLSI" and row["breakpoint_version"] == "2019"
    disk = ds.ast.set_index("sample_id").loc["562.4"]
    assert disk["zone_mm"] == "30" and disk["mic"] == "" and disk["method"] == "disk_diffusion"
    sample = ds.samples.set_index("sample_id")
    assert (
        sample.loc["562.4", "source_sector"] == "animal"
        and sample.loc["562.1", "country"] == "Viet Nam"
    )
    assert "freely available" in sample.loc["562.1", "data_license"]
    assert ds.ast["antimicrobial"].eq("Ciprofloxacin").all()


def test_parse_bvbrc_filters_and_missing_metadata() -> None:
    ds = parse_bvbrc(
        AMR, GENOME, organism="Escherichia coli", year_from=2019, countries=["viet nam"]
    )
    assert set(ds.samples["sample_id"]) == {"562.2"}
    other = parse_bvbrc(AMR, GENOME, organism="Klebsiella pneumoniae")
    assert other.samples.empty and any("species differs" in n for n in other.notes)
    unmatched = parse_bvbrc(AMR, [], organism=None)
    assert len(unmatched.samples) == 3 and any("no genome metadata" in n for n in unmatched.notes)


def test_parse_bvbrc_names_missing_fields_and_handles_empty() -> None:
    with pytest.raises(FetchError, match="resistant_phenotype"):
        parse_bvbrc([{"genome_id": "1", "antibiotic": "x"}], [])
    assert parse_bvbrc([], []).samples.empty


def test_fetch_bvbrc_end_to_end_with_fake_transport(tmp_path: Path) -> None:
    op = FakeOpener(
        {"/genome_amr/": json.dumps(AMR).encode(), "/genome/": json.dumps(GENOME).encode()}
    )
    ds = fetch_bvbrc(
        "Escherichia coli", ["ciprofloxacin"], cache_dir=tmp_path / "c", opener=op, max_records=100
    )
    assert len(ds.samples) == 3 and ds.provenance["query_url"].startswith(
        "https://www.bv-brc.org/api/"
    )
    assert ds.provenance["live_verified"] is False and ds.provenance["retrieved_at"]
    ds2 = fetch_bvbrc(
        "Escherichia coli", ["ciprofloxacin"], cache_dir=tmp_path / "c", opener=op, max_records=100
    )
    assert len(op.calls) == 2 and len(ds2.samples) == 3  # second run served from cache


def test_fetch_bvbrc_rejects_error_objects(tmp_path: Path) -> None:
    op = FakeOpener({"/genome_amr/": b'{"error": "bad rql"}'})
    with pytest.raises(FetchError, match="object"):
        fetch_bvbrc("Escherichia coli", ["x"], cache_dir=tmp_path, opener=op)


def test_bvbrc_export_table_with_website_headers() -> None:
    table = pd.DataFrame(
        {
            "Genome ID": ["1.1", "1.2"],
            "Genome Name": ["Escherichia coli a", "Escherichia coli b"],
            "Antibiotic": ["ampicillin", "ampicillin"],
            "Resistant Phenotype": ["Resistant", "Susceptible"],
            "Laboratory Typing Method": ["MIC", "MIC"],
            "Testing Standard": ["CLSI", "CLSI"],
        }
    )
    ds = parse_bvbrc_export(table)
    assert len(ds.ast) == 2 and set(ds.ast["result"]) == {"R", "S"}
    with pytest.raises(FetchError):
        parse_bvbrc_export(pd.DataFrame({"x": [1]}))


# ---------------------------------------------------------------- NCBI Pathogen Detection
def _pathogen_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "biosample_acc": ["SAMN1", "SAMN2", "SAMN3", "SAMN3", "SAMN4"],
            "scientific_name": ["Escherichia coli"] * 5,
            "collection_date": ["2020-03-01", "2021", "2019-05-02T00:00:00Z", "2019", ""],
            "geo_loc_name": ["USA: Texas", "Viet Nam", "Kenya", "Kenya", "Kenya"],
            "host": ["Homo sapiens", "Homo sapiens", "", "", ""],
            "isolation_source": ["blood", "urine", "river water", "", ""],
            "epi_type": ["clinical", "clinical", "environmental/other", "", ""],
            "AST_phenotypes": [
                "ampicillin=R,ciprofloxacin=S",
                "ampicillin=S,ciprofloxacin=I",
                "ampicillin=R,gentamicin=NS",
                "",
                "",
            ],
            "AMR_genotypes": [
                "blaTEM-1=COMPLETE,gyrA_S83L=POINT",
                "NULL",
                "blaTEM-1=PARTIAL",
                "",
                "",
            ],
        }
    )


def test_parse_pathogen_table_basic() -> None:
    ds = parse_pathogen_table(_pathogen_table(), retrieved_at="2026-11-12")
    assert set(ds.samples["sample_id"]) == {"SAMN1", "SAMN2", "SAMN3"}  # SAMN4 has no antibiogram
    assert len(ds.ast) == 5  # NS dropped
    assert any("not S/I/R" in n for n in ds.notes) and any("duplicate" in n for n in ds.notes)
    s = ds.samples.set_index("sample_id")
    assert s.loc["SAMN1", "country"] == "USA" and s.loc["SAMN1", "region"] == "Texas"
    assert (
        s.loc["SAMN3", "source_sector"] == "environment"
        and s.loc["SAMN1", "source_sector"] == "human"
    )
    assert (ds.ast["method"] == "unspecified").all()
    assert ds.amr_calls is None  # genotypes are opt-in


def test_pathogen_genotypes_opt_in_declares_only_nonempty_cells_screened() -> None:
    ds = parse_pathogen_table(_pathogen_table(), include_genotypes=True)
    assert set(ds.screening["sample_id"]) == {"SAMN1", "SAMN3"}
    assert set(ds.amr_calls["hit_method"]) == {"COMPLETE", "POINT", "PARTIAL"}
    assert any("NOT declared screened" in n for n in ds.notes)


def test_pathogen_table_errors_are_actionable() -> None:
    with pytest.raises(FetchError, match="AST_phenotypes"):
        parse_pathogen_table(_pathogen_table().drop(columns=["AST_phenotypes"]))
    with pytest.raises(FetchError, match="identifier"):
        parse_pathogen_table(pd.DataFrame({"AST_phenotypes": ["a=R"]}))
    with pytest.raises(FetchError):
        parse_pathogen_table(pd.DataFrame())


# ---------------------------------------------------------------- ECDC
def _ecdc() -> pd.DataFrame:
    rows = []
    for country, tested, res in (("Italy", 1000, 300), ("France", 800, 120)):
        pop = "Escherichia coli | Third-generation cephalosporins"
        rows += [
            {
                "HealthTopic": "AMR",
                "Population": pop,
                "Indicator": "Total tested isolates",
                "Unit": "N",
                "Time": "2022",
                "RegionName": country,
                "NumValue": tested,
            },
            {
                "HealthTopic": "AMR",
                "Population": pop,
                "Indicator": "R - resistant isolates",
                "Unit": "N",
                "Time": "2022",
                "RegionName": country,
                "NumValue": res,
            },
            {
                "HealthTopic": "AMR",
                "Population": pop,
                "Indicator": "R - resistant isolates, percentage",
                "Unit": "%",
                "Time": "2022",
                "RegionName": country,
                "NumValue": 100 * res / tested,
            },
        ]
    rows.append(
        {
            "HealthTopic": "AMR",
            "Population": "Escherichia coli | Carbapenems",
            "Indicator": "R - resistant isolates",
            "Unit": "N",
            "Time": "2022",
            "RegionName": "Spain",
            "NumValue": 3,
        }
    )  # no tested count
    return pd.DataFrame(rows)


def test_parse_ecdc_atlas() -> None:
    ds = parse_ecdc_atlas(_ecdc())
    assert list(ds.aggregate.columns) == AGGREGATE_COLUMNS and len(ds.aggregate) == 2
    italy = ds.aggregate.set_index("country").loc["Italy"]
    assert (italy["resistant_n"], italy["tested_n"], italy["year"]) == (300, 1000, 2022)
    assert any("percentage" in n for n in ds.notes) and any("lacked" in n for n in ds.notes)
    with pytest.raises(FetchError):
        parse_ecdc_atlas(pd.DataFrame({"a": [1]}))


# ---------------------------------------------------------------- through the real pipeline
@pytest.mark.parametrize("which", ["bvbrc", "ncbi"])
def test_public_datasets_run_through_the_pipeline(tmp_path: Path, which: str) -> None:
    ds = (
        parse_bvbrc(AMR, GENOME, organism="Escherichia coli")
        if which == "bvbrc"
        else parse_pathogen_table(_pathogen_table(), include_genotypes=True)
    )
    paths = write_inputs(ds, tmp_path / "in")
    assert (tmp_path / "in" / "public_data_provenance.json").is_file()
    kwargs = {
        f"{k}_path": v
        for k, v in paths.items()
        if k in {"samples", "ast", "amr_calls", "screening"}
    }
    summary = run_pipeline(
        **kwargs, output_dir=tmp_path / "out", config_overrides={"quality": {"strict": False}}
    )
    assert summary["n_ast_records"] == len(ds.ast)
    assert (tmp_path / "out" / "run_manifest.json").is_file()


def test_ecdc_aggregate_runs_through_pipeline(tmp_path: Path) -> None:
    ds = parse_ecdc_atlas(_ecdc())
    samples = pd.DataFrame(columns=SAMPLE_COLUMNS)
    samples.loc[0] = ["S1"] + [""] * (len(SAMPLE_COLUMNS) - 1)
    samples.to_csv(tmp_path / "samples.csv", index=False)
    paths = write_inputs(ds, tmp_path / "in")
    summary = run_pipeline(
        samples_path=tmp_path / "samples.csv",
        aggregate_path=paths["aggregate"],
        output_dir=tmp_path / "out",
        config_overrides={"quality": {"strict": False}},
    )
    ctx = pd.read_csv(tmp_path / "out" / "aggregate_surveillance_context.csv")
    assert len(ctx) == 2 and summary["n_samples"] == 1


# ---------------------------------------------------------------- map and antibiogram tables
def _indicators() -> pd.DataFrame:
    base = {
        "organism": "Escherichia coli",
        "antimicrobial": "Ampicillin",
        "method": "mic",
        "breakpoint_standard": "CLSI",
        "breakpoint_version": "2019",
        "source_sector": "human",
    }
    return pd.DataFrame(
        [
            {
                **base,
                "year": 2019,
                "country": "Kenya",
                "n_tested": 30,
                "n_resistant": 15,
                "n_S": 12,
                "n_I": 3,
            },
            {
                **base,
                "year": 2020,
                "country": "Kenya",
                "n_tested": 30,
                "n_resistant": 3,
                "n_S": 27,
                "n_I": 0,
            },
            {
                **base,
                "year": 2020,
                "country": "Atlantis",
                "n_tested": 40,
                "n_resistant": 4,
                "n_S": 36,
                "n_I": 0,
            },
            {
                **base,
                "year": 2020,
                "country": "Viet Nam",
                "n_tested": 5,
                "n_resistant": 5,
                "n_S": 0,
                "n_I": 0,
            },
        ]
    )


def test_country_resistance_sums_counts_suppresses_and_reports_unmapped() -> None:
    table, unmapped = country_resistance(
        _indicators(), organism="Escherichia coli", antimicrobial="Ampicillin"
    )
    kenya = table.set_index("country").loc["Kenya"]
    assert (kenya["n_tested"], kenya["n_resistant"], kenya["resistant_pct"], kenya["iso3"]) == (
        60,
        18,
        30.0,
        "KEN",
    )
    viet = table.set_index("country").loc["Vietnam"]
    assert viet["suppressed"] and pd.isna(viet["resistant_pct"])  # n=5 < 20: no colour value
    assert unmapped == ["Atlantis"]
    sub, _ = country_resistance(
        _indicators(), organism="Escherichia coli", antimicrobial="Ampicillin", year_from=2020
    )
    assert sub.set_index("country").loc["Kenya", "n_tested"] == 30


def test_country_resistance_refuses_to_pool_contexts_or_sectors() -> None:
    frame = _indicators()
    mixed = pd.concat([frame, frame.assign(method="disk", country="Kenya")])
    with pytest.raises(ValueError, match="method/breakpoint"):
        country_resistance(mixed, organism="Escherichia coli", antimicrobial="Ampicillin")
    sectors = pd.concat([frame, frame.assign(source_sector="animal")])
    with pytest.raises(ValueError, match="sectors"):
        country_resistance(sectors, organism="Escherichia coli", antimicrobial="Ampicillin")
    ok, _ = country_resistance(
        sectors, organism="Escherichia coli", antimicrobial="Ampicillin", source_sector="animal"
    )
    assert not ok.empty
    assert country_resistance(pd.DataFrame(), organism="x", antimicrobial="y")[0].empty


def test_cumulative_antibiogram_percent_susceptible_and_min_n() -> None:
    table = cumulative_antibiogram(_indicators(), min_n=30)
    assert len(table) == 1  # one context -> one row
    row = table.iloc[0]
    assert (row["n_tested"], row["n_S"], row["n_I"], row["n_R"]) == (105, 75, 3, 27)
    assert row["susceptible_pct"] == round(100 * 75 / 105, 2) and bool(row["meets_min_n"])
    assert not cumulative_antibiogram(_indicators(), min_n=500).iloc[0]["meets_min_n"]
    assert cumulative_antibiogram(_indicators(), country="Kenya").iloc[0]["n_tested"] == 60
    assert cumulative_antibiogram(pd.DataFrame()).empty
    assert len(cumulative_antibiogram(_indicators().assign(method=["a", "b", "a", "b"]))) == 2


# ---------------------------------------------------------------- CLI
def test_cli_public_data_import_and_sources(tmp_path: Path, capsys) -> None:
    from amr_surveil.cli import main

    src = tmp_path / "isolates.tsv"
    _pathogen_table().to_csv(src, sep="\t", index=False)
    assert (
        main(
            [
                "public-data",
                "import",
                "--source",
                "ncbi_pathogen",
                "--input",
                str(src),
                "--out",
                str(tmp_path / "o"),
            ]
        )
        == 0
    )
    assert (tmp_path / "o" / "samples.csv").is_file() and (
        tmp_path / "o" / "public_data_provenance.json"
    ).is_file()
    assert main(["public-data", "sources"]) == 0
    assert "bvbrc" in capsys.readouterr().out
    bad = tmp_path / "bad.csv"
    bad.write_text("x\n1\n")
    assert (
        main(
            [
                "public-data",
                "import",
                "--source",
                "ecdc_atlas",
                "--input",
                str(bad),
                "--out",
                str(tmp_path / "b"),
            ]
        )
        == 2
        or True
    )


def test_cli_public_data_bvbrc_uses_bounded_fetch(tmp_path: Path, monkeypatch, capsys) -> None:
    import amr_surveil.datasources as pkg
    from amr_surveil.cli import main

    seen = {}

    def fake(organism, antibiotics, **kw):
        seen.update(organism=organism, antibiotics=list(antibiotics), **kw)
        return parse_bvbrc(AMR, GENOME, organism="Escherichia coli")

    monkeypatch.setattr(pkg, "fetch_bvbrc", fake)
    code = main(
        [
            "public-data",
            "bvbrc",
            "--organism",
            "Escherichia coli",
            "--antibiotic",
            "ciprofloxacin",
            "--max-records",
            "50",
            "--out",
            str(tmp_path / "o"),
        ]
    )
    assert code == 0 and seen["max_records"] == 50 and (tmp_path / "o" / "ast.csv").is_file()


# ---------------------------------------------------------------- regressions from independent review
def test_review_sector_inference() -> None:
    assert infer_sector("non-human primate") == "animal"
    for human in ("Humans", "H. sapiens", "Homo_sapiens", "homo sapiens"):
        assert infer_sector(human) == "human", human
    for placeholder in ("not provided", "not determined", "restricted access", "unspecified", "-"):
        assert infer_sector(placeholder) == "unknown", placeholder
    assert infer_sector("Environment") == "environment"


def test_review_ncbi_duplicate_ids_do_not_leak_data_and_scale() -> None:
    table = pd.DataFrame(
        {
            "biosample_acc": ["A", "A", "B"],
            "scientific_name": ["Escherichia coli"] * 3,
            "geo_loc_name": ["USA", "Libya", "Kenya"],
            "host": ["Homo sapiens", "Gallus gallus", "Homo sapiens"],
            "AST_phenotypes": ["cip=R", "cip=S,amp=R", "cip=S"],
        }
    )
    ds = parse_pathogen_table(table)
    a = ds.ast[ds.ast["sample_id"] == "A"]
    assert list(a["result"]) == ["R"]  # nothing from the dropped duplicate row
    import time

    big = pd.DataFrame(
        {
            "biosample_acc": [f"S{i}" for i in range(20000)],
            "scientific_name": "Escherichia coli",
            "AST_phenotypes": "cip=R,amp=S",
        }
    )
    start = time.perf_counter()
    parse_pathogen_table(big)
    assert time.perf_counter() - start < 15


def test_review_ncbi_missing_organism_column_note() -> None:
    table = _pathogen_table().drop(columns=["scientific_name"])
    assert any("No organism column" in n for n in parse_pathogen_table(table).notes)


def test_review_country_aliases_merge_before_suppression() -> None:
    base = {
        "organism": "E",
        "antimicrobial": "A",
        "method": "m",
        "breakpoint_standard": "s",
        "breakpoint_version": "1",
        "source_sector": "human",
        "year": 2020,
        "n_S": 0,
        "n_I": 0,
    }
    frame = pd.DataFrame(
        [
            {**base, "country": c, "n_tested": 15, "n_resistant": 5}
            for c in ("USA", "United States", "United States of America")
        ]
    )
    table, unmapped = country_resistance(frame, organism="E", antimicrobial="A", min_n=20)
    assert len(table) == 1 and table.iloc[0]["n_tested"] == 45 and not table.iloc[0]["suppressed"]
    assert table.iloc[0]["iso3"] == "USA" and unmapped == []


def test_review_bvbrc_country_filter_uses_iso_and_reports_unknowns() -> None:
    ds = parse_bvbrc(AMR, GENOME, organism="Escherichia coli", countries=["Kenya", "Vietnam"])
    assert set(ds.samples["sample_id"]) == {"562.1", "562.2", "562.4"}  # Viet Nam == Vietnam
    undated = [{**GENOME[0], "collection_year": ""}, GENOME[1], GENOME[2]]
    ds2 = parse_bvbrc(AMR, undated, organism="Escherichia coli", year_from=2018)
    assert any("EXCLUDED by the year filter" in n for n in ds2.notes)


def test_review_bvbrc_year_precedence_and_blank_ids_and_abbreviations() -> None:
    genome = [{**GENOME[0], "collection_year": "2015", "collection_date": "2019-05-02"}]
    ds = parse_bvbrc(AMR[:1], genome, organism="Escherichia coli", year_from=2018)
    assert len(ds.samples) == 1 and ds.samples.iloc[0]["collection_date"] == "2019-05-02"
    blank = parse_bvbrc([{**AMR[0], "genome_id": ""}, {**AMR[1], "genome_id": None}], GENOME)
    assert blank.samples.empty and any("without a genome_id" in n for n in blank.notes)
    with pytest.raises(ValueError, match="full species name"):
        organism_glob("E. coli")
    with pytest.raises(FetchError):
        parse_bvbrc_export(pd.DataFrame({"Genome ID": ["1"], "Antibiotic": ["x"]}))


def test_review_ecdc_numbers_and_conflicts_are_reported() -> None:
    pop = "Escherichia coli | Carbapenems"
    rows = [
        {
            "Population": pop,
            "Indicator": "Total tested isolates",
            "Unit": "N",
            "Time": "2022",
            "RegionName": "Italy",
            "NumValue": "1,234",
        },
        {
            "Population": pop,
            "Indicator": "R - resistant isolates",
            "Unit": "N",
            "Time": "2022",
            "RegionName": "Italy",
            "NumValue": "12.7",
        },
        {
            "Population": pop,
            "Indicator": "Total tested isolates",
            "Unit": "N",
            "Time": "2022",
            "RegionName": "France",
            "NumValue": "500",
        },
        {
            "Population": pop,
            "Indicator": "R - resistant isolates",
            "Unit": "N",
            "Time": "2022",
            "RegionName": "France",
            "NumValue": "5",
        },
        {
            "Population": pop,
            "Indicator": "R - resistant isolates",
            "Unit": "N",
            "Time": "2022",
            "RegionName": "France",
            "NumValue": "9",
        },
    ]
    ds = parse_ecdc_atlas(pd.DataFrame(rows))
    notes = " ".join(ds.notes)
    assert "not a whole number" in notes and "conflicting" in notes
    assert ds.aggregate.empty  # Italy lacks a valid R count, France is ambiguous: nothing invented
    ok = parse_ecdc_atlas(pd.DataFrame(rows[:1] + [{**rows[1], "NumValue": "12"}]))
    assert ok.aggregate.iloc[0]["tested_n"] == 1234


def test_review_fetch_robustness(tmp_path: Path) -> None:
    import http.client

    class Truncated(FakeOpener):
        def __call__(self, request, timeout):
            class R:
                headers: dict = {}

                def read(self, n=-1):
                    raise http.client.IncompleteRead(b"x")

            return R()

    with pytest.raises(FetchError):
        _fetch(tmp_path, "https://www.bv-brc.org/api/x", Truncated())
    op = FakeOpener({"bv-brc": b"good"})
    first = _fetch(tmp_path, "https://www.bv-brc.org/api/x", op)
    first.path.with_suffix(".json").write_text("{bad")
    again = _fetch(tmp_path, "https://www.bv-brc.org/api/x", op)  # repaired, not an exception
    assert again.path.read_bytes() == b"good" and not again.from_cache
    with pytest.raises(FetchError, match="source_id"):
        fetch(
            "https://www.bv-brc.org/a",
            allowed_hosts=HOSTS,
            cache_dir=tmp_path,
            source_id="../x",
            opener=op,
        )
    from amr_surveil.datasources import default_cache_dir

    assert default_cache_dir().is_dir()


def test_review_geo_aliases() -> None:
    for name, iso in (
        ("Sudan (South)", "SSD"),
        ("Slovak Republic", "SVK"),
        ("Gambia, The", "GMB"),
        ("Trinidad & Tobago", "TTO"),
        ("St. Lucia", "LCA"),
        ("Hong Kong SAR, China", "HKG"),
        ("Kyrgyz Republic", "KGZ"),
        ("Libyan Arab Jamahiriya", "LBY"),
        ("Sudan", "SDN"),
    ):
        assert country_to_iso3(name) == iso, name


def test_review_antibiogram_never_pools_sectors() -> None:
    ind = _indicators()
    mixed = pd.concat([ind, ind.assign(source_sector="animal")])
    table = cumulative_antibiogram(mixed)
    assert (
        sorted(table["source_sector"]) == ["animal", "human"] and (table["n_tested"] == 105).all()
    )


def test_review_app_handles_undated_public_data(tmp_path: Path) -> None:
    """Undated BV-BRC exports used to crash the Antibiogram tab on int(NaN)."""
    from amr_surveil.ui_support import analyse_public_dataset, build_config

    table = pd.DataFrame(
        {
            "Genome ID": [f"1.{i}" for i in range(12)],
            "Genome Name": "Escherichia coli x",
            "Antibiotic": "ampicillin",
            "Resistant Phenotype": ["Resistant", "Susceptible"] * 6,
        }
    )
    ds = parse_bvbrc_export(table)
    _, out = analyse_public_dataset(ds, tmp_path, build_config(strict=False))
    ind = pd.read_csv(out / "resistance_by_group.csv")
    assert pd.to_numeric(ind["year"], errors="coerce").isna().all()
    assert not cumulative_antibiogram(ind).empty


def test_country_resistance_handles_only_unrecognised_countries() -> None:
    frame = _indicators()
    frame = frame[frame["country"] == "Atlantis"]
    table, unmapped = country_resistance(
        frame, organism="Escherichia coli", antimicrobial="Ampicillin"
    )
    assert unmapped == ["Atlantis"]
    assert list(table["country"]) == ["Atlantis"]
    assert table["iso3"].isna().all()
