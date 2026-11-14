"""FishAMR-Link adapters.

Three ways to bring FishAMR-Link evidence into AMR-Surveil, in order of preference:

1. **AMRFinderPlus reports** produced by (or alongside) FishAMR-Link: use
   :mod:`amr_surveil.adapters.amrfinder`. This is the validated detection engine.
2. :func:`normalize_fishamr_feature_table` for the in-app "Genomics Core" table
   (columns ``gene, class, mechanism, feature_type, confidence``), with an explicit evidence tier
   inferred from how each row was produced.
3. :func:`read_fishamr_warehouse` for the SQLite warehouse (``isolates`` /
   ``genomic_sequences``), excluding the seeded reference baseline by default.
4. :func:`normalize_fishamr_calls` / :func:`normalize_fishamr_plasmids`: generic JSON column
   mapping for any other CSV export.

Findings from inspecting FishAMR-Link v17.0 (commit fa84114) that shape this adapter
------------------------------------------------------------------------------------
* The app has no CSV export of its own. Rows carry no sample identifier, so one must be supplied.
* Its fallback scanner matches short peptide motifs and reports fixed ``confidence`` constants
  (0.89 / 0.92); native AMRFinderPlus rows carry 1.0. These are not statistical confidences, so
  the value is ignored and an ``evidence_tier`` is recorded instead.
* The remote-BLAST fallback records generic nt homology hits labelled ``Homology Match``.
* The warehouse's seeded baseline uses simulated sequences, and ``phenotypic_resistance`` is a
  free-text drug list with no S/I/R category, method or breakpoint. It is therefore **never**
  converted into AST records.
* Latitude/longitude are intentionally not exported (avoid exposing precise
  locations in public repositories).
"""

from __future__ import annotations

import json
import re
import sqlite3
from pathlib import Path

import pandas as pd

from amr_surveil.adapters.common import map_columns
from amr_surveil.constants import AMR_CALL_COLUMNS, PLASMID_COLUMNS, SAMPLE_COLUMNS
from amr_surveil.exceptions import ConfigurationError, InputFileError
from amr_surveil.io import read_csv, write_csv

AMR_TARGETS = set(AMR_CALL_COLUMNS)
PLASMID_TARGETS = set(PLASMID_COLUMNS)

_SEEDED_BASELINE = {
    "NZ_CP012345.1",
    "NZ_CP098765.1",
    "NZ_CP112233.1",
    "NZ_CP445566.1",
    "NZ_CP778899.1",
    "NZ_CP881122.1",
    "NZ_CP334455.1",
    "NZ_CP990011.1",
    "NZ_CP556677.1",
    "NZ_CP223344.1",
}
SEEDED_BASELINE_ACCESSIONS = _SEEDED_BASELINE | {a + "_v" for a in _SEEDED_BASELINE}

_AMRFINDER_METHOD = re.compile(
    r"^(ALLELE|EXACT|BLAST|PARTIAL|HMM|POINT|INTERNAL_STOP|COMPLETE)[A-Z_]*$"
)
_PARTIAL_METHOD = re.compile(r"^(PARTIAL|INTERNAL_STOP)")
_FALLBACK_MECHANISMS = {
    "Enzymatic Inactivation",
    "Active Efflux Sluice System",
    "Dihydropteroate Synthase Bypass",
    "LPS Modification Cluster",
    "Enzymatic Modification",
    "Replicative Transposition",
    "Site-Specific Recombination",
    "Autonomous Replication Initiation",
    "Iron Acquisition Sluice",
    "Pore-Forming Cytolytic Activity",
}


def infer_fishamr_tier(gene: object, klass: object, mechanism: object) -> str:
    """Evidence tier for one in-app FishAMR-Link feature row (see module docstring)."""
    mech = "" if pd.isna(mechanism) else str(mechanism).strip()
    cls = "" if pd.isna(klass) else str(klass).strip()
    if mech.startswith("Identified via remote blastn") or cls == "Homology Match":
        return "homology_search"
    if _PARTIAL_METHOD.match(mech):
        return "homology_search"  # truncated / incomplete hit: not a functional determinant
    if _AMRFINDER_METHOD.match(mech):
        return "curated_database"
    if mech in _FALLBACK_MECHANISMS:
        return "heuristic_motif"
    return "unspecified"


def normalize_fishamr_calls(
    input_path: str | Path, output_path: str | Path, mapping_path: str | Path
) -> pd.DataFrame:
    """Normalize any FishAMR-Link-like gene-call CSV using a JSON column mapping."""
    return map_columns(
        input_path,
        output_path,
        mapping_path,
        required_targets={"sample_id", "gene"},
        allowed_targets=AMR_TARGETS,
        output_columns=AMR_CALL_COLUMNS,
        default_source_system="FishAMR-Link",
    )


def normalize_fishamr_plasmids(
    input_path: str | Path, output_path: str | Path, mapping_path: str | Path
) -> pd.DataFrame:
    """Normalize FishAMR-Link plasmid/replicon exports using a JSON column mapping."""
    return map_columns(
        input_path,
        output_path,
        mapping_path,
        required_targets={"sample_id", "replicon"},
        allowed_targets=PLASMID_TARGETS,
        output_columns=PLASMID_COLUMNS,
        default_source_system="FishAMR-Link",
    )


def normalize_fishamr_feature_table(
    input_path: str | Path,
    output_calls: str | Path,
    output_plasmids: str | Path | None = None,
    *,
    sample_id: str | None = None,
    organism: str | None = None,
) -> dict[str, object]:
    """Convert the in-app ``amr_genes`` table (CSV/TSV) into canonical calls and replicon rows."""
    path = Path(input_path)
    frame = read_csv(
        path,
        table="fishamr_feature_table",
        sep="\t" if path.suffix.lower() in {".tsv", ".txt"} else ",",
    )
    for needed in ("gene", "feature_type"):
        if needed not in frame.columns:
            raise InputFileError(
                f"{path}: missing column '{needed}' (expected the FishAMR-Link feature table)"
            )
    if "sample_id" not in frame.columns:
        if not sample_id:
            raise ConfigurationError(
                "The FishAMR-Link feature table has no sample_id column; pass sample_id explicitly"
            )
        frame["sample_id"] = sample_id
    for col in ("class", "mechanism"):
        if col not in frame.columns:
            frame[col] = pd.NA
    frame["evidence_tier"] = [
        infer_fishamr_tier(g, c, m)
        for g, c, m in zip(frame["gene"], frame["class"], frame["mechanism"], strict=True)
    ]
    ftype = frame["feature_type"].astype("string").str.strip()
    ftype_key = ftype.str.casefold()
    is_amr = (ftype_key == "amr_gene").fillna(False).astype(bool)
    is_replicon = (ftype_key == "mge").fillna(False).astype(bool) & frame["gene"].astype(
        "string"
    ).str.contains(r"(?:^Inc|rep[A-Z]?$|_rep)", regex=True, na=False)
    method = (
        frame["evidence_tier"]
        .map(
            {
                "curated_database": "AMRFinderPlus (via FishAMR-Link)",
                "homology_search": "NCBI remote BLASTn (via FishAMR-Link)",
                "heuristic_motif": "FishAMR-Link peptide-motif scan",
            }
        )
        .fillna("FishAMR-Link")
    )
    calls = pd.DataFrame(
        {
            "sample_id": frame["sample_id"],
            "organism": organism,
            "gene": frame["gene"],
            "method": method,
            "source_system": "FishAMR-Link",
            "element_type": ftype,
            "determinant_class": frame["class"],
            "evidence_tier": frame["evidence_tier"],
        }
    )[is_amr.to_numpy()]
    for col in AMR_CALL_COLUMNS:
        if col not in calls.columns:
            calls[col] = pd.NA
    calls = calls[AMR_CALL_COLUMNS]
    write_csv(calls, output_calls)
    plasmids = pd.DataFrame(
        {
            "sample_id": frame["sample_id"],
            "replicon": frame["gene"],
            "mobility_prediction": pd.NA,
            "method": method,
            "database": pd.NA,
            "database_version": pd.NA,
            "source_system": "FishAMR-Link",
            "evidence_tier": frame["evidence_tier"],
        }
    )[is_replicon.to_numpy()][PLASMID_COLUMNS]
    if output_plasmids is not None:
        write_csv(plasmids, output_plasmids)
    other = ftype[~is_amr & ~is_replicon].fillna("missing").value_counts().to_dict()
    return {
        "n_rows": len(frame),
        "n_amr_calls": len(calls),
        "n_replicon_rows": len(plasmids),
        "excluded_feature_types": {str(k): int(v) for k, v in other.items()},
        "calls_by_evidence_tier": calls["evidence_tier"].value_counts().to_dict(),
    }


def read_fishamr_warehouse(
    db_path: str | Path,
    samples_out: str | Path,
    calls_out: str | Path,
    *,
    include_seeded_baseline: bool = False,
) -> dict[str, object]:
    """Export isolates and detected genes from a FishAMR-Link SQLite warehouse (read-only)."""
    path = Path(db_path)
    if not path.is_file():
        raise InputFileError(f"Warehouse not found: {path}")
    try:
        con = sqlite3.connect(f"{path.resolve().as_uri()}?mode=ro", uri=True)
    except sqlite3.Error as exc:
        raise InputFileError(f"Cannot open warehouse {path}: {exc}") from exc
    try:
        tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if not {"isolates", "genomic_sequences"} <= tables:
            raise InputFileError(
                f"{path}: expected tables 'isolates' and 'genomic_sequences', found {sorted(tables)}"
            )
        rows = con.execute(
            "SELECT i.accession, i.organism, i.host_species, i.country, i.collection_year, i.isolation_source, "
            "g.detected_amr_genes FROM isolates i LEFT JOIN genomic_sequences g ON g.isolate_id = i.id "
            "ORDER BY i.accession, g.id"
        ).fetchall()
    except sqlite3.Error as exc:
        raise InputFileError(f"Cannot read warehouse {path}: {exc}") from exc
    finally:
        con.close()
    sample_rows, call_rows, seen = [], [], set()
    n_seeded = 0
    for accession, organism, host, country, year, source, genes_json in rows:
        seeded = accession in SEEDED_BASELINE_ACCESSIONS
        if seeded and not include_seeded_baseline:
            n_seeded += 1 if accession not in seen else 0
            seen.add(accession)
            continue
        if accession not in seen:
            seen.add(accession)
            sample_rows.append(
                {
                    "sample_id": accession,
                    "organism": organism,
                    "collection_date": None if year is None else str(year),
                    "country": country,
                    "source_sector": "aquaculture",
                    "specimen": source,
                    "host": host,
                    "sequencing_accession": accession,
                    "source_study": "FishAMR-Link seeded baseline (simulated sequences)"
                    if seeded
                    else "FishAMR-Link warehouse",
                    "data_license": "synthetic" if seeded else None,
                }
            )
        for gene in json.loads(genes_json) if genes_json else []:
            symbol = gene.get("gene") if isinstance(gene, dict) else gene
            if symbol:
                call_rows.append(
                    {
                        "sample_id": accession,
                        "organism": organism,
                        "gene": str(symbol),
                        "method": "FishAMR-Link warehouse record",
                        "source_system": "FishAMR-Link",
                        "evidence_tier": "unspecified",
                    }
                )
    samples = pd.DataFrame(sample_rows, columns=SAMPLE_COLUMNS)
    calls = pd.DataFrame(call_rows, columns=AMR_CALL_COLUMNS)
    write_csv(samples, samples_out)
    write_csv(calls, calls_out)
    return {
        "n_isolates_exported": len(samples),
        "n_calls_exported": len(calls),
        "n_seeded_baseline_isolates_skipped": n_seeded,
        "included_seeded_baseline": include_seeded_baseline,
        "note": "Warehouse phenotypic_resistance is free text without S/I/R, method or breakpoint and is not exported as AST.",
    }
