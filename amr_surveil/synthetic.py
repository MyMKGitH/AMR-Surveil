"""Deterministic, clearly synthetic data for demos, tests and benchmarks.

Nothing here is real surveillance data. The generator deliberately plants structure so every
analysis path has something to find, and documents it so results can be checked against ground
truth:

* **Planted genotype-phenotype association.** ``blaCTX-M-15``->Cefotaxime, ``qnrS1``->Ciprofloxacin,
  ``blaKPC-2``->Meropenem. A carried gene gives resistance with probability 0.92, an absent gene
  with probability 0.06 (``PLANTED``). ``tetA``/``sul1`` are decoys with no mapped drug.
* **Method change in 2025.** *E. coli* cefotaxime in 2025 is tested by disk diffusion (zone only),
  earlier years by broth microdilution (MIC) - a comparability threat the report must flag.
* **Repeat patients** (same patient, species, year) so first-isolate deduplication has an effect.
* **Low-evidence calls.** A few aquatic samples carry ``heuristic_motif`` calls that the default
  configuration excludes.
* Synthetic breakpoints are labelled ``SYNTHETIC`` and are **not** EUCAST or CLSI values.

With ``inject_defects=True`` a documented set of data-quality defects is added so the audit and
the strict quality gate can be exercised (see :data:`DEFECTS`).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

PLANTED = {"p_resistant_if_gene": 0.92, "p_resistant_if_no_gene": 0.06}
GENE_FOR_DRUG = {"Cefotaxime": "blaCTX-M-15", "Ciprofloxacin": "qnrS1", "Meropenem": "blaKPC-2"}
DRUGS = list(GENE_FOR_DRUG)
SYNTHETIC_BREAKPOINTS = {  # drug: (mic_s_le, mic_r_gt, zone_s_ge, zone_r_lt)
    "Cefotaxime": (1.0, 2.0, 20.0, 17.0),
    "Ciprofloxacin": (0.25, 0.5, 25.0, 22.0),
    "Meropenem": (2.0, 8.0, 22.0, 16.0),
}
MIC_LADDER = [0.03, 0.06, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
DEFECTS = [
    "duplicate_sample_row",
    "non_iso_collection_date",
    "future_collection_date",
    "organism_case_variant",
    "invalid_ast_category",
    "blank_ast_result",
    "orphan_ast_sample_id",
    "duplicate_aggregate_key",
    "low_identity_call",
    "mic_unparseable",
]


def _fmt(x: float) -> str:
    return f"{x:g}"


def generate_synthetic_dataset(
    output_dir: str | Path,
    n_samples: int = 120,
    seed: int = 20261107,
    *,
    inject_defects: bool = False,
) -> dict[str, Path]:
    """Write synthetic input tables and return their paths (keys match the pipeline inputs)."""
    if n_samples < 8:
        raise ValueError("n_samples must be at least 8")
    rng = np.random.default_rng(seed)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    ids = [f"SYN-{i:05d}" for i in range(1, n_samples + 1)]
    organism = np.where(np.arange(n_samples) % 3 == 0, "Escherichia coli", "Klebsiella pneumoniae")
    sector = np.where(np.arange(n_samples) % 4 == 0, "aquatic", "human")
    years = rng.choice([2022, 2023, 2024, 2025], size=n_samples)
    months = rng.integers(1, 13, n_samples)
    days = rng.integers(1, 28, n_samples)
    patient = [f"PAT-{i:05d}" if s == "human" else None for i, s in enumerate(sector, start=1)]
    dates = [f"{y}-{int(m):02d}-{int(d):02d}" for y, m, d in zip(years, months, days, strict=True)]
    for i in range(1, n_samples):  # repeat patients: same patient/species/year, later date
        if sector[i] == "human" and sector[i - 1] == "human" and i % 9 == 5:
            patient[i] = patient[i - 1]
            organism[i] = organism[i - 1]
            years[i] = years[i - 1]
            dates[i] = f"{years[i]}-12-{int(rng.integers(10, 28)):02d}"
    specimen = [
        ("blood" if i % 2 else "urine")
        if s == "human"
        else ("surface_water" if i % 2 else "fish_tissue")
        for i, s in enumerate(sector)
    ]
    samples = pd.DataFrame(
        {
            "sample_id": ids,
            "organism": organism,
            "collection_date": dates,
            "country": "Exampleland",
            "region": np.where(sector == "human", "Central", "River Basin"),
            "source_sector": sector,
            "specimen": specimen,
            "patient_id": patient,
            "host": np.where(sector == "human", "human", "environmental_sample"),
            "sequencing_accession": [f"SYNACC{i:06d}" for i in range(1, n_samples + 1)],
            "sequencing_type": "Illumina short-read (synthetic)",
            "source_study": "AMR-Surveil synthetic fixture",
            "data_license": "CC0-1.0 synthetic data",
        }
    )

    # ---- genotype: which of the three mapped genes + two decoys each sample carries
    carried: dict[str, set[str]] = {sid: set() for sid in ids}
    for sid in ids:
        for gene in [*GENE_FOR_DRUG.values(), "tetA", "sul1"]:
            if rng.random() < (0.30 if gene in GENE_FOR_DRUG.values() else 0.25):
                carried[sid].add(gene)
    call_rows = []
    for i, sid in enumerate(ids):
        for gene in sorted(carried[sid]):
            tier = "curated_database"
            method = "AMRFinderPlus"
            if sector[i] == "aquatic" and gene in {"tetA", "sul1"} and rng.random() < 0.5:
                tier, method = "heuristic_motif", "FishAMR-Link peptide-motif scan"
            call_rows.append(
                {
                    "sample_id": sid,
                    "organism": organism[i],
                    "gene": gene,
                    "identity_pct": _fmt(round(float(rng.uniform(97.0, 100.0)), 2)),
                    "coverage_pct": _fmt(round(float(rng.uniform(90.0, 100.0)), 2)),
                    "contig": f"contig_{int(rng.integers(1, 20))}",
                    "start": "1000",
                    "end": "1800",
                    "method": method,
                    "database": "SYNTHETIC-AMR-DB",
                    "database_version": "2026.1",
                    "source_system": "FishAMR-Link" if sector[i] == "aquatic" else "AMRFinderPlus",
                    "element_type": "AMR",
                    "determinant_class": "SYNTHETIC",
                    "hit_method": "EXACTX" if tier == "curated_database" else "MOTIF",
                    "evidence_tier": tier,
                }
            )
    amr_calls = pd.DataFrame(call_rows)
    screening = pd.DataFrame(
        {
            "sample_id": ids,
            "method": "AMRFinderPlus (synthetic)",
            "tool_version": "synthetic-3.12",
            "database": "SYNTHETIC-AMR-DB",
            "database_version": "2026.1",
            "identity_threshold_pct": "90",
            "coverage_threshold_pct": "50",
            "source_system": [
                "FishAMR-Link" if s == "aquatic" else "AMRFinderPlus" for s in sector
            ],
        }
    )

    # ---- phenotype from MIC via the synthetic breakpoints (so recorded categories are consistent)
    ast_rows = []
    for i, sid in enumerate(ids):
        for drug in DRUGS:
            s_le, r_gt, z_s, z_r = SYNTHETIC_BREAKPOINTS[drug]
            has_gene = GENE_FOR_DRUG[drug] in carried[sid]
            resistant = rng.random() < (
                PLANTED["p_resistant_if_gene"] if has_gene else PLANTED["p_resistant_if_no_gene"]
            )
            if resistant:
                mic = float(rng.choice([m for m in MIC_LADDER if m > r_gt]))
            else:
                mic = float(rng.choice([m for m in MIC_LADDER if m <= r_gt]))
            category = "S" if mic <= s_le else ("R" if mic > r_gt else "I")
            row = {
                "sample_id": sid,
                "organism": organism[i],
                "antimicrobial": drug,
                "result": category,
                "mic": _fmt(mic),
                "mic_unit": "mg/L",
                "zone_mm": None,
                "method": "broth_microdilution",
                "breakpoint_standard": "SYNTHETIC",
                "breakpoint_version": "0.0",
                "interpretation_date": "2026-01-01",
                "testing_lab": "Synthetic Lab A" if i % 2 == 0 else "Synthetic Lab B",
            }
            if years[i] == 2025 and organism[i] == "Escherichia coli" and drug == "Cefotaxime":
                zone = (
                    float(rng.integers(10, 17)) if category == "R" else float(rng.integers(22, 31))
                )
                if category == "I":
                    zone = 18.0
                row.update(
                    mic=None,
                    mic_unit=None,
                    zone_mm=_fmt(zone),
                    method="disk_diffusion",
                    result="S" if zone >= z_s else ("R" if zone < z_r else "I"),
                )
            ast_rows.append(row)
    ast = pd.DataFrame(ast_rows)

    breakpoints = pd.DataFrame(
        [
            {
                "organism": org,
                "antimicrobial": drug,
                "standard": "SYNTHETIC",
                "version": "0.0",
                "mic_s_le": _fmt(v[0]),
                "mic_r_gt": _fmt(v[1]),
                "zone_s_ge": _fmt(v[2]),
                "zone_r_lt": _fmt(v[3]),
                "source_url": "https://example.invalid/synthetic-breakpoints",
            }
            for org in ("Escherichia coli", "Klebsiella pneumoniae")
            for drug, v in SYNTHETIC_BREAKPOINTS.items()
        ]
    )
    gene_map = pd.DataFrame(
        [
            {
                "organism": org,
                "gene": gene,
                "antimicrobial": drug,
                "mapping_source": "Synthetic fixture mapping (not a curated database)",
                "match_type": "exact",
            }
            for org in ("Escherichia coli", "Klebsiella pneumoniae")
            for drug, gene in GENE_FOR_DRUG.items()
        ]
    )

    aggregate = pd.DataFrame(
        [
            {
                "year": y,
                "country": "Exampleland",
                "source_sector": sec,
                "organism": "Escherichia coli",
                "antimicrobial": "Cefotaxime",
                "resistant_n": r,
                "tested_n": t,
                "source_name": "Synthetic aggregate fixture",
                "metric_definition": "R isolates / isolates tested for this organism-antimicrobial group",
                "source_url": "https://example.invalid/synthetic",
            }
            for y, sec, r, t in (
                ("2023", "human", "18", "100"),
                ("2024", "human", "22", "110"),
                ("2024", "aquatic", "6", "20"),
            )
        ]
    )
    incidence = pd.DataFrame(
        [
            {
                "year": "2024",
                "country": "Exampleland",
                "source_sector": "human",
                "organism": "Escherichia coli",
                "antimicrobial": "Cefotaxime",
                "resistant_infection_n": "8",
                "population_n": "1000000",
                "population_multiplier": "100000",
                "source_name": "Synthetic incidence fixture",
                "case_definition": "Synthetic resistant bloodstream infection episodes",
                "source_url": "https://example.invalid/synthetic-incidence",
            }
        ]
    )
    consumption = pd.DataFrame(
        [
            {
                "year": "2024",
                "country": "Exampleland",
                "care_sector": "community",
                "antimicrobial": "Amoxicillin",
                "atc_code": "J01CA04",
                "ddd_total": "1000000",
                "population_n": "100000",
                "period_days": "365",
                "aware_group": "Access",
                "source_name": "Synthetic consumption fixture",
                "data_source_type": "dispensing",
                "source_url": "https://example.invalid/synthetic-consumption",
            }
        ]
    )

    pred_rows = []
    for i, sid in enumerate(ids[: n_samples // 2]):
        drug = DRUGS[i % len(DRUGS)]
        gene_present = GENE_FOR_DRUG[drug] in carried[sid]
        predicted_r = gene_present if rng.random() < 0.85 else not gene_present
        pred_rows.append(
            {
                "sample_id": sid,
                "organism": organism[i],
                "antimicrobial": drug,
                "predicted_mic": None,
                "predicted_mic_unit": None,
                "predicted_class": "predicted_resistant"
                if predicted_r
                else "predicted_susceptible",
                "prediction_score": _fmt(round(float(rng.uniform(0.6, 0.99)), 3)),
                "model_name": "EpiResNet",
                "model_version": "synthetic-example",
                "checkpoint_id": "not-a-real-checkpoint",
                "training_data_version": "synthetic-only",
                "split_strategy": "synthetic-example-no-performance-claim",
                "source_system": "EpiResNet",
                "prediction_unit": "isolate-antibiotic pair",
                "protein_cluster": None,
                "decision_threshold": "0.5",
                "probability_temperature": "1.0",
                "evaluation_split": "test",
            }
        )
    predictions = pd.DataFrame(pred_rows)
    plasmids = pd.DataFrame(
        [
            {
                "sample_id": ids[i],
                "replicon": "IncFIB-like",
                "mobility_prediction": "unknown",
                "method": "synthetic_placeholder",
                "database": "SYNTHETIC-PLASMID-DB",
                "database_version": "0.0",
                "source_system": "AMRFinderPlus",
                "evidence_tier": "curated_database",
            }
            for i in range(0, n_samples, max(1, n_samples // 8))
        ]
    )

    if inject_defects:
        samples = pd.concat([samples, samples.iloc[[0]]], ignore_index=True)  # duplicate_sample_row
        samples.loc[1, "collection_date"] = "03/05/2024"  # non_iso_collection_date
        samples.loc[2, "collection_date"] = "2099-01-01"  # future_collection_date
        samples.loc[3, "organism"] = "escherichia coli "  # organism_case_variant
        ast.loc[0, "result"] = "X"  # invalid_ast_category
        ast.loc[1, "result"] = None  # blank_ast_result
        ast = pd.concat(
            [ast, ast.iloc[[2]].assign(sample_id="NOT-IN-MANIFEST")], ignore_index=True
        )  # orphan
        aggregate = pd.concat(
            [aggregate, aggregate.iloc[[0]]], ignore_index=True
        )  # duplicate_aggregate_key
        amr_calls.loc[0, "identity_pct"] = "62.5"  # low_identity_call
        ast.loc[3, "mic"] = "1,5"  # mic_unparseable

    tables = {
        "samples": ("samples.csv", samples),
        "amr_calls": ("amr_calls.csv", amr_calls),
        "screening": ("screening.csv", screening),
        "ast": ("ast.csv", ast),
        "breakpoints": ("breakpoints.csv", breakpoints),
        "gene_drug_map": ("gene_drug_map.csv", gene_map),
        "aggregate": ("aggregate_surveillance.csv", aggregate),
        "predictions": ("model_predictions.csv", predictions),
        "plasmids": ("plasmid_calls.csv", plasmids),
        "incidence": ("incidence.csv", incidence),
        "consumption": ("consumption.csv", consumption),
    }
    paths: dict[str, Path] = {}
    for key, (name, frame) in tables.items():
        path = out / name
        frame.to_csv(path, index=False)
        paths[key] = path
    return paths
