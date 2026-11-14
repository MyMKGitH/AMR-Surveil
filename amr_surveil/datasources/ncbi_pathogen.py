"""NCBI Pathogen Detection isolate tables (Isolates Browser TSV / metadata TSV) -> AMR-Surveil inputs.

The user downloads the file (NCBI's own terms say it cannot grant rights on behalf of submitters, so
this project never bundles or redistributes it). Columns are matched loosely by header name and the
mapping that was used is returned so the user can check it.

``AST_phenotypes`` is read as ``drug=R,drug=S,...`` (BioSample antibiogram calls). Method and
breakpoint standard are normally absent from this table and are recorded as ``unspecified``.
"""

from __future__ import annotations

import re

import pandas as pd

from amr_surveil.constants import AMR_CALL_COLUMNS, AST_COLUMNS, SAMPLE_COLUMNS, SCREENING_COLUMNS
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
from amr_surveil.datasources.fetch import FetchError
from amr_surveil.datasources.registry import get_source, license_text

ID_ALIASES = (
    "biosample_acc",
    "BioSample",
    "biosample",
    "#label",
    "label",
    "asm_acc",
    "Assembly",
    "target_acc",
)
ORGANISM_ALIASES = ("scientific_name", "Scientific name", "Organism group", "taxgroup_name")
DATE_ALIASES = ("collection_date", "Collection date")
COUNTRY_ALIASES = ("geo_loc_name", "Location", "Geographic location")
HOST_ALIASES = ("host", "Host")
SOURCE_ALIASES = ("isolation_source", "Isolation source")
EPI_ALIASES = ("epi_type", "Isolation type")
AST_ALIASES = ("AST_phenotypes", "AST phenotypes")
GENO_ALIASES = ("AMR_genotypes", "AMR genotypes")
ASSEMBLY_ALIASES = ("asm_acc", "Assembly", "biosample_acc", "BioSample")

_PAIR = re.compile(r"^\s*(?P<key>[^=]+?)\s*=\s*(?P<val>[^=]+?)\s*$")


def parse_pairs(text: object) -> list[tuple[str, str]]:
    """``a=R,b=S`` (or ``;``-separated) -> [(a, R), (b, S)]; malformed fragments are skipped."""
    value = clean(text)
    if not value:
        return []
    out = []
    for part in re.split(r"[;,]", value):
        match = _PAIR.match(part)
        if match:
            out.append((match.group("key"), match.group("val")))
    return out


def parse_pathogen_table(
    table: pd.DataFrame, *, include_genotypes: bool = False, retrieved_at: str = ""
) -> PublicDataset:
    if table.empty:
        raise FetchError("The table has no rows.")
    table = table.reset_index(drop=True)
    ids = pick(table, *ID_ALIASES)
    ast_col = pick(table, *AST_ALIASES)
    if (ids.map(clean) == "").all():
        raise FetchError(
            "No isolate identifier column found (looked for biosample_acc / BioSample / #label / asm_acc). "
            f"Columns present: {list(table.columns)[:10]}."
        )
    present = {header_key(c) for c in table.columns}
    if not any(header_key(a) in present for a in AST_ALIASES):
        raise FetchError(
            "The table has no AST_phenotypes column, so there is nothing to analyse as susceptibility data. "
            "In the Isolates Browser, keep the 'AST phenotypes' column and filter to isolates that have it."
        )
    notes: list[str] = []
    organism = pick(table, *ORGANISM_ALIASES).map(clean)
    if not any(header_key(a) in present for a in ORGANISM_ALIASES):
        notes.append("No organism column found; organism is blank and will be flagged by QA.")
    elif not any(header_key(a) in present for a in ORGANISM_ALIASES[:2]):
        notes.append(
            "Only an organism-group column was found; organism is the group name, not the species."
        )
    host = pick(table, *HOST_ALIASES).map(clean)
    source = pick(table, *SOURCE_ALIASES).map(clean)
    epi = pick(table, *EPI_ALIASES).map(clean)
    frame = pd.DataFrame(
        {
            "sample_id": ids.map(clean),
            "organism": organism,
            "collection_date": pick(table, *DATE_ALIASES).map(normalise_date),
            "country": pick(table, *COUNTRY_ALIASES).map(
                lambda v: clean(v).split(":", 1)[0].strip()
            ),
            "region": pick(table, *COUNTRY_ALIASES).map(
                lambda v: clean(v).split(":", 1)[1].strip() if ":" in clean(v) else ""
            ),
            "source_sector": [
                infer_sector(h, s, e) for h, s, e in zip(host, source, epi, strict=True)
            ],
            "specimen": source,
            "patient_id": "",
            "host": host,
            "sequencing_accession": pick(table, *ASSEMBLY_ALIASES).map(clean),
            "sequencing_type": "",
            "source_study": "NCBI Pathogen Detection (user-downloaded table)",
            "data_license": license_text("ncbi_pathogen", retrieved_at or "date not recorded"),
        }
    )
    blank_id = frame["sample_id"] == ""
    if blank_id.any():
        notes.append(f"{int(blank_id.sum())} row(s) without an identifier were dropped.")
    frame = frame[~blank_id]
    dup = frame["sample_id"].duplicated(keep="first")
    if dup.any():
        notes.append(f"{int(dup.sum())} duplicate identifier row(s) dropped (first kept).")
    frame = frame[~dup]
    first_rows = list(frame.index)  # only the kept (first) row of each isolate may contribute data

    ast_rows: list[dict[str, str]] = []
    skipped_pairs = 0
    skipped_values: dict[str, int] = {}
    for idx in first_rows:
        sid, text = frame.at[idx, "sample_id"], ast_col.at[idx]
        for drug, value in parse_pairs(text):
            result = map_result(value)
            if result is None:
                skipped_pairs += 1
                skipped_values[value] = skipped_values.get(value, 0) + 1
                continue
            ast_rows.append({"sample_id": sid, "antimicrobial": title_drug(drug), "result": result})
    if skipped_pairs:
        notes.append(
            f"{skipped_pairs} drug result(s) skipped (not S/I/R): {dict(list(skipped_values.items())[:6])}."
        )
    org_by_id = dict(zip(frame["sample_id"], frame["organism"], strict=True))
    ast = pd.DataFrame(ast_rows, columns=["sample_id", "antimicrobial", "result"])
    ast["organism"] = ast["sample_id"].map(org_by_id)
    for col, default in (
        ("mic", ""),
        ("mic_unit", ""),
        ("zone_mm", ""),
        ("method", "unspecified"),
        ("breakpoint_standard", "unspecified"),
        ("breakpoint_version", ""),
        ("interpretation_date", ""),
        ("testing_lab", ""),
    ):
        ast[col] = default
    ast = ast[AST_COLUMNS].drop_duplicates()
    with_ast = set(ast["sample_id"])
    frame = frame[frame["sample_id"].isin(with_ast)]
    notes.append(
        f"{len(frame)} isolate(s) carry AST results; isolates without an antibiogram were left out."
    )
    notes.append(
        "Method and breakpoint standard are not in this table (recorded as 'unspecified'); results from "
        "different laboratories and standards are therefore of unknown comparability."
    )
    notes.append(
        "Pathogen Detection is a convenience sample (outbreak-enriched); no patient ID exists, so "
        "deduplication is impossible. Proportions are not prevalence."
    )
    dataset = PublicDataset(
        "ncbi_pathogen",
        frame[SAMPLE_COLUMNS].reset_index(drop=True),
        ast.reset_index(drop=True),
        {
            "source_id": "ncbi_pathogen",
            "source_name": get_source("ncbi_pathogen").name,
            "terms": get_source("ncbi_pathogen").terms,
            "terms_url": get_source("ncbi_pathogen").terms_url,
            "terms_checked": get_source("ncbi_pathogen").terms_checked,
            "attribution": get_source("ncbi_pathogen").attribution,
            "retrieved_at": retrieved_at,
            "live_verified": False,
            "include_genotypes": include_genotypes,
        },
        notes,
    )
    if include_genotypes:
        _add_genotypes(dataset, table, ids, first_rows)
    return dataset


def _add_genotypes(
    dataset: PublicDataset, table: pd.DataFrame, ids: pd.Series, first_rows: list[int]
) -> None:
    geno_col = pick(table, *GENO_ALIASES)
    keep = set(dataset.samples["sample_id"])
    org_by_id = dict(zip(dataset.samples["sample_id"], dataset.samples["organism"], strict=True))
    rows, screened = [], set()
    for idx in first_rows:
        sid, text = clean(ids.at[idx]), geno_col.at[idx]
        if sid not in keep or not clean(text):
            continue
        screened.add(sid)
        for gene, status in parse_pairs(text):
            rows.append(
                {
                    "sample_id": sid,
                    "organism": org_by_id.get(sid, ""),
                    "gene": gene,
                    "identity_pct": "",
                    "coverage_pct": "",
                    "contig": "",
                    "start": "",
                    "end": "",
                    "method": "AMRFinderPlus (NCBI Pathogen Detection)",
                    "database": "unspecified",
                    "database_version": "unspecified",
                    "source_system": "NCBI Pathogen Detection",
                    "element_type": "AMR",
                    "determinant_class": "",
                    "hit_method": status.upper(),
                    "evidence_tier": "curated_database",
                }
            )
    dataset.amr_calls = pd.DataFrame(rows, columns=AMR_CALL_COLUMNS)
    dataset.screening = pd.DataFrame(
        [
            {
                "sample_id": sid,
                "method": "AMRFinderPlus (NCBI Pathogen Detection)",
                "tool_version": "unspecified",
                "database": "unspecified",
                "database_version": "unspecified",
                "identity_threshold_pct": "",
                "coverage_threshold_pct": "",
                "source_system": "NCBI Pathogen Detection",
            }
            for sid in sorted(screened)
        ],
        columns=SCREENING_COLUMNS,
    )
    dataset.notes.append(
        f"Genotypes imported for {len(screened)} isolate(s) with a non-empty AMR_genotypes cell. Isolates "
        "with an empty cell were NOT declared screened (empty can mean 'not run'), which biases "
        "genotype-negative counts downward — treat concordance as exploratory."
    )
