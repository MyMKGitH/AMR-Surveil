"""BV-BRC (Bacterial and Viral Bioinformatics Resource Center) AMR phenotype connector.

Live path: query the public Data API (``genome_amr`` for lab AST records, then ``genome`` for country,
year, host), join them, and convert to canonical AMR-Surveil inputs. File path: convert a
"AMR phenotypes" table exported from the BV-BRC website (headers are matched loosely).

The API base URL and RQL syntax follow the published BV-BRC API documentation; **field names could
not be tested against the live server by the maintainers**, so every missing expected field is
reported by name instead of failing silently (see :func:`parse_bvbrc`).
"""

from __future__ import annotations

import json
import urllib.parse
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd

from amr_surveil.constants import AST_COLUMNS, SAMPLE_COLUMNS
from amr_surveil.datasources.common import (
    PublicDataset,
    clean,
    header_key,
    infer_sector,
    map_result,
    normalise_date,
    pick,
    title_drug,
)
from amr_surveil.datasources.fetch import FetchError, Opener, fetch
from amr_surveil.datasources.geo import country_to_iso3
from amr_surveil.datasources.registry import get_source, license_text

API_BASE = "https://www.bv-brc.org/api"
AMR_FIELDS = (
    "genome_id,genome_name,antibiotic,resistant_phenotype,measurement,measurement_sign,"
    "measurement_value,measurement_unit,laboratory_typing_method,testing_standard,testing_standard_year"
)
GENOME_FIELDS = (
    "genome_id,species,isolation_country,collection_year,collection_date,host_name,"
    "isolation_source,biosample_accession,assembly_accession"
)
ID_CHUNK = 100
MAX_RECORDS_CEILING = 25_000


def _enc(value: str) -> str:
    return urllib.parse.quote(str(value), safe="*")


def organism_glob(organism: str) -> str:
    """Wildcard pattern for the genome name, e.g. 'Escherichia coli' -> 'Escherichia*coli*'."""
    parts = [p for p in organism.replace(".", " ").split() if p]
    if not parts:
        raise ValueError("organism must not be empty")
    if "." in organism.split()[0] or len(parts[0]) < 3:
        raise ValueError(
            f"Use the full species name (e.g. 'Escherichia coli'), not an abbreviation like {organism!r}."
        )
    return "*".join(_enc(p) for p in parts) + "*"


def amr_query_url(organism: str, antibiotics: Iterable[str], *, limit: int, offset: int = 0) -> str:
    drugs = [d.strip().lower() for d in antibiotics if d and d.strip()]
    if not drugs:
        raise ValueError("choose at least one antibiotic")
    if not 1 <= limit <= MAX_RECORDS_CEILING:
        raise ValueError(f"limit must be between 1 and {MAX_RECORDS_CEILING}")
    drug_clause = (
        f"eq(antibiotic,{_enc(drugs[0])})"
        if len(drugs) == 1
        else "in(antibiotic,(" + ",".join(_enc(d) for d in drugs) + "))"
    )
    name_clause = f"eq(genome_name,{organism_glob(organism)})"
    return (
        f"{API_BASE}/genome_amr/?and({drug_clause},{name_clause})"
        f"&select({AMR_FIELDS})&limit({limit},{offset})"
    )


def genome_query_url(genome_ids: Iterable[str]) -> str:
    ids = ",".join(_enc(g) for g in genome_ids)
    return f"{API_BASE}/genome/?in(genome_id,({ids}))&select({GENOME_FIELDS})&limit({ID_CHUNK})"


def _load_records(path: Path) -> list[dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(data, dict):
        raise FetchError(
            f"BV-BRC returned an object instead of a list ({str(data)[:200]}). The query may be invalid."
        )
    if not isinstance(data, list):
        raise FetchError("BV-BRC returned an unexpected response shape.")
    return [r for r in data if isinstance(r, dict)]


def _method(typing: object, platform: object = "") -> str:
    text = clean(typing).lower().replace(" ", "_")
    return text or "unspecified"


def parse_bvbrc(
    amr_records: list[dict[str, Any]],
    genome_records: list[dict[str, Any]],
    *,
    organism: str | None = None,
    year_from: int | None = None,
    year_to: int | None = None,
    countries: Iterable[str] | None = None,
    retrieved_at: str = "",
) -> PublicDataset:
    """Convert BV-BRC ``genome_amr`` + ``genome`` records to canonical samples/AST tables."""
    notes: list[str] = []
    amr = pd.DataFrame(amr_records)
    genome = pd.DataFrame(genome_records)
    expected_amr = {"genome_id", "antibiotic", "resistant_phenotype"}
    missing = sorted(expected_amr - set(amr.columns)) if not amr.empty else []
    if missing:
        raise FetchError(
            f"BV-BRC AMR records lack expected field(s) {missing}; got {sorted(amr.columns)[:12]}. "
            "The API schema may have changed — please open an issue with this message."
        )
    if amr.empty:
        return PublicDataset(
            "bvbrc",
            pd.DataFrame(columns=SAMPLE_COLUMNS),
            pd.DataFrame(columns=AST_COLUMNS),
            _provenance(retrieved_at, organism, {}),
            ["The query returned no AMR records."],
        )
    n_raw = len(amr)
    blank_id = amr["genome_id"].map(lambda v: clean(v) == "")
    if blank_id.any():
        notes.append(f"{int(blank_id.sum())} AST record(s) without a genome_id were dropped.")
        amr = amr.loc[~blank_id].copy()
    amr["genome_id"] = amr["genome_id"].astype(str)
    mapped = amr["resistant_phenotype"].map(map_result)
    dropped = amr.loc[mapped.isna(), "resistant_phenotype"].map(lambda v: clean(v) or "(blank)")
    if len(dropped):
        counts = dropped.value_counts().head(6).to_dict()
        notes.append(
            f"{int(len(dropped))} of {n_raw} AST record(s) dropped because the phenotype is not S/I/R: {counts}. "
            "'Non-susceptible' is not recoded because it merges I and R."
        )
    amr = amr.loc[mapped.notna()].copy()
    amr["result"] = mapped.dropna()

    meta = genome.copy() if not genome.empty else pd.DataFrame(columns=["genome_id"])
    if "genome_id" in meta.columns:
        meta["genome_id"] = meta["genome_id"].astype(str)
        meta = meta.drop_duplicates("genome_id")
    for column in (
        "species",
        "isolation_country",
        "collection_year",
        "collection_date",
        "host_name",
        "isolation_source",
        "biosample_accession",
        "assembly_accession",
    ):
        if column not in meta.columns:
            meta[column] = ""
    no_meta = ~amr["genome_id"].isin(meta["genome_id"])
    if no_meta.any():
        notes.append(
            f"{int(no_meta.sum())} AST record(s) have no genome metadata (country/year unknown); kept with blanks."
        )
    merged = amr.merge(meta, on="genome_id", how="left")

    if organism:
        wanted = organism.strip().lower()
        species = merged["species"].map(lambda v: clean(v).lower())
        keep = (species == wanted) | (species == "")
        if (species != "").any():
            removed = int((~keep).sum())
            if removed:
                notes.append(f"{removed} record(s) removed: species differs from '{organism}'.")
            merged = merged[keep]
    # Same precedence as the sample's collection_date: a full date wins over the year field.
    year = pd.to_numeric(merged["collection_date"].map(normalise_date).str[:4], errors="coerce")
    year = year.fillna(
        pd.to_numeric(merged["collection_year"].map(normalise_date), errors="coerce")
    )
    if year_from is not None or year_to is not None:
        keep_year = pd.Series(True, index=merged.index)
        if year_from is not None:
            keep_year &= year >= year_from
        if year_to is not None:
            keep_year &= year <= year_to
        undated = int(year.isna().sum())
        outside = int((~keep_year & year.notna()).sum())
        if undated:
            notes.append(
                f"{undated} AST record(s) have no collection year and were EXCLUDED by the year filter."
            )
        if outside:
            notes.append(
                f"{outside} AST record(s) outside {year_from or '…'}-{year_to or '…'} removed."
            )
        merged = merged[keep_year.fillna(False)]
    if countries:
        wanted_iso = {country_to_iso3(c) for c in countries if c.strip()} - {None}
        wanted_names = {c.strip().lower() for c in countries if c.strip()}

        def _wanted(value: object) -> bool:
            name = clean(value)
            iso = country_to_iso3(name)
            return (iso is not None and iso in wanted_iso) or name.lower() in wanted_names

        keep_country = merged["isolation_country"].map(_wanted)
        n_no_country = int(merged["isolation_country"].map(lambda v: clean(v) == "").sum())
        if n_no_country:
            notes.append(
                f"{n_no_country} AST record(s) have no country and were EXCLUDED by the country filter."
            )
        merged = merged[keep_country]
    merged = merged.reset_index(drop=True)

    sample_rows = merged.drop_duplicates("genome_id")
    org_fallback = (
        sample_rows["genome_name"].map(lambda v: " ".join(clean(v).split()[:2]))
        if "genome_name" in sample_rows.columns
        else ""
    )
    species_col = sample_rows["species"].map(clean)
    samples = pd.DataFrame(
        {
            "sample_id": sample_rows["genome_id"],
            "organism": species_col.where(species_col != "", org_fallback),
            "collection_date": sample_rows.apply(
                lambda r: (
                    normalise_date(r["collection_date"]) or normalise_date(r["collection_year"])
                ),
                axis=1,
            ),
            "country": sample_rows["isolation_country"].map(clean),
            "region": "",
            "source_sector": sample_rows.apply(
                lambda r: infer_sector(r["host_name"], r["isolation_source"]), axis=1
            ),
            "specimen": sample_rows["isolation_source"].map(clean),
            "patient_id": "",
            "host": sample_rows["host_name"].map(clean),
            "sequencing_accession": sample_rows["assembly_accession"]
            .map(clean)
            .where(lambda s: s != "", sample_rows["biosample_accession"].map(clean)),
            "sequencing_type": "",
            "source_study": "BV-BRC genome_amr (public AMR phenotypes)",
            "data_license": license_text("bvbrc", retrieved_at or "date not recorded"),
        }
    )[SAMPLE_COLUMNS]

    org_by_id = dict(zip(samples["sample_id"], samples["organism"], strict=True))
    typing = merged.get("laboratory_typing_method", pd.Series([""] * len(merged)))
    method = typing.map(_method)
    is_zone = method.str.contains("disk|disc|diffusion", regex=True)
    value = merged.get("measurement_value", pd.Series([""] * len(merged))).map(clean)
    sign = merged.get("measurement_sign", pd.Series([""] * len(merged))).map(clean)
    unit = merged.get("measurement_unit", pd.Series([""] * len(merged))).map(clean)
    standard = merged.get("testing_standard", pd.Series([""] * len(merged))).map(clean)
    std_year = merged.get("testing_standard_year", pd.Series([""] * len(merged))).map(clean)
    ast = pd.DataFrame(
        {
            "sample_id": merged["genome_id"],
            "organism": merged["genome_id"].map(org_by_id),
            "antimicrobial": merged["antibiotic"].map(title_drug),
            "result": merged["result"],
            "mic": (sign + value).where(~is_zone & (value != ""), ""),
            "mic_unit": unit.where(~is_zone & (value != ""), ""),
            "zone_mm": value.where(is_zone, ""),
            "method": method,
            "breakpoint_standard": standard.where(standard != "", "unspecified"),
            "breakpoint_version": std_year,
            "interpretation_date": "",
            "testing_lab": "",
        }
    )[AST_COLUMNS]
    ast = ast.drop_duplicates()
    notes.append(
        "Public-repository genomes are a convenience sample (outbreak and resistant isolates are "
        "over-represented); no patient identifier exists, so isolates cannot be de-duplicated per patient."
    )
    return PublicDataset(
        "bvbrc",
        samples,
        ast,
        _provenance(retrieved_at, organism, {"n_raw_amr_records": n_raw}),
        notes,
    )


def _provenance(retrieved_at: str, organism: str | None, extra: dict[str, Any]) -> dict[str, Any]:
    info = get_source("bvbrc")
    return {
        "source_id": "bvbrc",
        "source_name": info.name,
        "terms": info.terms,
        "terms_url": info.terms_url,
        "terms_checked": info.terms_checked,
        "attribution": info.attribution,
        "retrieved_at": retrieved_at,
        "organism_query": organism,
        "live_verified": info.live_verified,
        **extra,
    }


def fetch_bvbrc(
    organism: str,
    antibiotics: Iterable[str],
    *,
    cache_dir: str | Path,
    max_records: int = 5000,
    year_from: int | None = None,
    year_to: int | None = None,
    countries: Iterable[str] | None = None,
    opener: Opener | None = None,
    refresh: bool = False,
) -> PublicDataset:
    """Query BV-BRC for lab AST records of ``organism`` and the listed antibiotics (bounded)."""
    hosts = get_source("bvbrc").hosts
    headers = {"Accept": "application/json"}
    url = amr_query_url(organism, antibiotics, limit=max_records)
    first = fetch(
        url,
        allowed_hosts=hosts,
        cache_dir=cache_dir,
        source_id="bvbrc_amr",
        headers=headers,
        opener=opener,
        refresh=refresh,
    )
    amr_records = _load_records(first.path)
    ids = sorted({str(r.get("genome_id")) for r in amr_records if r.get("genome_id")})
    genome_records: list[dict[str, Any]] = []
    stamps = [first.retrieved_at]
    digests = {"genome_amr": first.sha256}
    for start in range(0, len(ids), ID_CHUNK):
        chunk = ids[start : start + ID_CHUNK]
        result = fetch(
            genome_query_url(chunk),
            allowed_hosts=hosts,
            cache_dir=cache_dir,
            source_id="bvbrc_genome",
            headers=headers,
            opener=opener,
            refresh=refresh,
        )
        genome_records.extend(_load_records(result.path))
        stamps.append(result.retrieved_at)
    dataset = parse_bvbrc(
        amr_records,
        genome_records,
        organism=organism,
        year_from=year_from,
        year_to=year_to,
        countries=countries,
        retrieved_at=max(stamps),
    )
    dataset.provenance.update(
        {
            "query_url": url,
            "response_sha256": digests,
            "n_genome_metadata_requests": -(-len(ids) // ID_CHUNK),
        }
    )
    if len(amr_records) >= max_records:
        dataset.notes.append(
            f"The result was cut at {max_records} records; narrow the query for a complete picture."
        )
    return dataset


def parse_bvbrc_export(table: pd.DataFrame, *, retrieved_at: str = "") -> PublicDataset:
    """Convert a BV-BRC 'AMR phenotypes' table downloaded from the website (CSV/TSV)."""
    key = {header_key(c) for c in table.columns}
    if not {"genomeid", "antibiotic", "resistantphenotype"} <= key:
        raise FetchError(
            "This does not look like a BV-BRC AMR phenotype export (needs 'Genome ID', 'Antibiotic' and "
            f"'Resistant Phenotype' columns; found {list(table.columns)[:8]})."
        )
    amr_records = pd.DataFrame(
        {
            "genome_id": pick(table, "genome_id", "Genome ID"),
            "genome_name": pick(table, "genome_name", "Genome Name"),
            "antibiotic": pick(table, "antibiotic"),
            "resistant_phenotype": pick(table, "resistant_phenotype", "Resistant Phenotype"),
            "measurement_sign": pick(table, "measurement_sign", "Measurement Sign"),
            "measurement_value": pick(table, "measurement_value", "Measurement Value"),
            "measurement_unit": pick(table, "measurement_unit", "Measurement Unit"),
            "laboratory_typing_method": pick(
                table, "laboratory_typing_method", "Laboratory Typing Method"
            ),
            "testing_standard": pick(table, "testing_standard", "Testing Standard"),
            "testing_standard_year": pick(table, "testing_standard_year", "Testing Standard Year"),
        }
    ).to_dict("records")
    genome_records = pd.DataFrame(
        {
            "genome_id": pick(table, "genome_id", "Genome ID"),
            "species": pick(table, "species", "Species"),
            "isolation_country": pick(table, "isolation_country", "Isolation Country"),
            "collection_year": pick(table, "collection_year", "Collection Year"),
            "collection_date": pick(table, "collection_date", "Collection Date"),
            "host_name": pick(table, "host_name", "Host Name"),
            "isolation_source": pick(table, "isolation_source", "Isolation Source"),
        }
    ).to_dict("records")
    dataset = parse_bvbrc(amr_records, genome_records, retrieved_at=retrieved_at)
    dataset.notes.append(
        "Imported from a file you downloaded. If the export has no country/year columns, the map and "
        "time trends will be empty for these records."
    )
    return dataset
