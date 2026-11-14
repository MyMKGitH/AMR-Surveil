# AMR-Surveil analysis report

**Software version:** `1.1.1`  
**Run timestamp (UTC):** `2026-11-14T12:19:59.140822+00:00`  
**Reference date for date checks:** `2026-11-07`  
**Deduplication rule:** `none`  
**Quality gate:** `strict`

> **Scope:** research and portfolio analysis of supplied data. Not a clinical decision-support system, a diagnostic report, a national surveillance submission, or evidence of operational validation.

## Study definition

- **Name:** AMR-Surveil synthetic demonstration
- **Question:** Does the pipeline recover the planted gene-phenotype association and flag the planted method change?
- **Inclusion criteria:** All synthetic isolates.
- **Exclusion criteria:** None.
- **Sample-selection description:** Synthetic data generated with a fixed seed; no real isolates or sampling frame exist.

## Dataset inventory

| Stream | Rows | Interpretation |
|---|---:|---|
| Isolate/sample manifest | 120 | Isolate-level metadata |
| Screened samples (genotyping declared) | 120 | Only these can have a *negative* genotype |
| AMR determinant calls (supplied / used) | 169 / 164 | Positive calls; evidence tiers excluded: `heuristic_motif, homology_search` |
| Phenotypic AST records | 360 | Observed categories supplied by the user |
| Aggregate surveillance rows | 3 | Source counts; not isolate-linked |
| Population incidence rows | 1 | Cases per stated denominator; not isolate-linked |
| Antimicrobial consumption rows | 1 | DDD-based indicator; not patient-level use |
| Model predictions | 60 | Research outputs; not observed AST |
| Plasmid/replicon rows | 8 | Upstream evidence; not proof of transfer |
| Isolates with AMR calls and AST (exact ID) | 87 | Co-occurrence inventory only |
| Records excluded with a recorded reason | 5 | See `exclusion_log.csv` |

## Data-quality review

Errors **0** · Warnings **1** · Info **0**

| Severity | Table | Code | Rows | Finding |
|---|---|---|---|---|
| WARNING | amr_calls | low_confidence_evidence_tier | 5 | Calls from heuristic-motif or generic homology-search methods are not validated resistance-gene detections. They are excluded from linkage and concordance by the configured evidence-tier policy. |

### Exclusion summary

| Stage | Table | Reason | Records |
|---|---|---|---|
| evidence_tier | amr_calls | excluded_evidence_tier | 5 |

## Isolate-level AST resistance indicators

Resistance proportion is **R / (S + I + R)** per group. `I` is kept in the denominator and not in the numerator (EUCAST: *susceptible, increased exposure*). Wilson score intervals. The small-n flag is a display warning, not a protocol rule.

| Year | Country | Sector | Organism | Antimicrobial | Breakpoint | R / tested | Resistance (95% CI) | Caution |
|---|---|---|---|---|---|---|---|---|
| 2025 | Exampleland | human | Klebsiella pneumoniae | Ciprofloxacin | SYNTHETIC 0.0 | 5/20 | 25.0% (11.2–46.9) | — |
| 2025 | Exampleland | human | Klebsiella pneumoniae | Meropenem | SYNTHETIC 0.0 | 7/20 | 35.0% (18.1–56.7) | — |
| 2025 | Exampleland | human | Klebsiella pneumoniae | Cefotaxime | SYNTHETIC 0.0 | 4/20 | 20.0% (8.1–41.6) | — |
| 2022 | Exampleland | human | Klebsiella pneumoniae | Ciprofloxacin | SYNTHETIC 0.0 | 10/18 | 55.6% (33.7–75.4) | — |
| 2022 | Exampleland | human | Klebsiella pneumoniae | Cefotaxime | SYNTHETIC 0.0 | 4/18 | 22.2% (9.0–45.2) | — |
| 2022 | Exampleland | human | Klebsiella pneumoniae | Meropenem | SYNTHETIC 0.0 | 7/18 | 38.9% (20.3–61.4) | — |
| 2023 | Exampleland | human | Klebsiella pneumoniae | Ciprofloxacin | SYNTHETIC 0.0 | 4/14 | 28.6% (11.7–54.6) | — |
| 2023 | Exampleland | human | Klebsiella pneumoniae | Cefotaxime | SYNTHETIC 0.0 | 3/14 | 21.4% (7.6–47.6) | — |
| 2023 | Exampleland | human | Klebsiella pneumoniae | Meropenem | SYNTHETIC 0.0 | 6/14 | 42.9% (21.4–67.4) | — |
| 2024 | Exampleland | aquatic | Klebsiella pneumoniae | Ciprofloxacin | SYNTHETIC 0.0 | 5/9 | 55.6% (26.7–81.1) | Small n |
| 2024 | Exampleland | aquatic | Klebsiella pneumoniae | Cefotaxime | SYNTHETIC 0.0 | 4/9 | 44.4% (18.9–73.3) | Small n |
| 2022 | Exampleland | human | Escherichia coli | Cefotaxime | SYNTHETIC 0.0 | 3/9 | 33.3% (12.1–64.6) | Small n |
| 2022 | Exampleland | human | Escherichia coli | Meropenem | SYNTHETIC 0.0 | 5/9 | 55.6% (26.7–81.1) | Small n |
| 2022 | Exampleland | human | Escherichia coli | Ciprofloxacin | SYNTHETIC 0.0 | 2/9 | 22.2% (6.3–54.7) | Small n |
| 2024 | Exampleland | aquatic | Klebsiella pneumoniae | Meropenem | SYNTHETIC 0.0 | 2/9 | 22.2% (6.3–54.7) | Small n |
| 2025 | Exampleland | human | Escherichia coli | Ciprofloxacin | SYNTHETIC 0.0 | 3/9 | 33.3% (12.1–64.6) | Small n |
| 2025 | Exampleland | human | Escherichia coli | Cefotaxime | SYNTHETIC 0.0 | 3/9 | 33.3% (12.1–64.6) | Small n |
| 2025 | Exampleland | human | Escherichia coli | Meropenem | SYNTHETIC 0.0 | 1/9 | 11.1% (2.0–43.5) | Small n |
| 2024 | Exampleland | human | Klebsiella pneumoniae | Ciprofloxacin | SYNTHETIC 0.0 | 1/8 | 12.5% (2.2–47.1) | Small n |
| 2024 | Exampleland | human | Klebsiella pneumoniae | Cefotaxime | SYNTHETIC 0.0 | 1/8 | 12.5% (2.2–47.1) | Small n |
| 2024 | Exampleland | human | Klebsiella pneumoniae | Meropenem | SYNTHETIC 0.0 | 3/8 | 37.5% (13.7–69.4) | Small n |
| 2024 | Exampleland | human | Escherichia coli | Meropenem | SYNTHETIC 0.0 | 1/6 | 16.7% (3.0–56.4) | Small n |
| 2023 | Exampleland | aquatic | Klebsiella pneumoniae | Ciprofloxacin | SYNTHETIC 0.0 | 2/6 | 33.3% (9.7–70.0) | Small n |
| 2023 | Exampleland | human | Escherichia coli | Ciprofloxacin | SYNTHETIC 0.0 | 2/6 | 33.3% (9.7–70.0) | Small n |
| 2024 | Exampleland | human | Escherichia coli | Cefotaxime | SYNTHETIC 0.0 | 1/6 | 16.7% (3.0–56.4) | Small n |

*Showing 25 of 48 rows; see the CSV output for all rows.*

## Comparability flags

Flags mean *do not read a trend into these numbers until the cause is understood*. Nothing is adjusted or smoothed.

| Severity | Organism | Antimicrobial | Country | Sector | Flag | Detail |
|---|---|---|---|---|---|---|
| WARNING | Escherichia coli | Cefotaxime | Exampleland | aquatic | small_denominator_years | years with fewer than 10 isolates: 2022 (n=3), 2023 (n=2), 2024 (n=2), 2025 (n=3) |
| WARNING | Escherichia coli | Cefotaxime | Exampleland | aquatic | ast_method_changes_across_years | AST method differs between years: 2022: broth_microdilution; 2023: broth_microdilution; 2024: broth_microdilution; 2025: disk_diffusion |
| WARNING | Escherichia coli | Ciprofloxacin | Exampleland | aquatic | small_denominator_years | years with fewer than 10 isolates: 2022 (n=3), 2023 (n=2), 2024 (n=2), 2025 (n=3) |
| WARNING | Escherichia coli | Meropenem | Exampleland | aquatic | small_denominator_years | years with fewer than 10 isolates: 2022 (n=3), 2023 (n=2), 2024 (n=2), 2025 (n=3) |
| WARNING | Klebsiella pneumoniae | Cefotaxime | Exampleland | aquatic | small_denominator_years | years with fewer than 10 isolates: 2022 (n=3), 2023 (n=6), 2024 (n=9), 2025 (n=2) |
| WARNING | Klebsiella pneumoniae | Cefotaxime | Exampleland | aquatic | testing_volume_jump | tested isolates changed from 3 (2022) to 6 (2023); possible change in panel, laboratory participation or sampling |
| WARNING | Klebsiella pneumoniae | Cefotaxime | Exampleland | aquatic | testing_volume_jump | tested isolates changed from 9 (2024) to 2 (2025); possible change in panel, laboratory participation or sampling |
| WARNING | Klebsiella pneumoniae | Ciprofloxacin | Exampleland | aquatic | small_denominator_years | years with fewer than 10 isolates: 2022 (n=3), 2023 (n=6), 2024 (n=9), 2025 (n=2) |
| WARNING | Klebsiella pneumoniae | Ciprofloxacin | Exampleland | aquatic | testing_volume_jump | tested isolates changed from 3 (2022) to 6 (2023); possible change in panel, laboratory participation or sampling |
| WARNING | Klebsiella pneumoniae | Ciprofloxacin | Exampleland | aquatic | testing_volume_jump | tested isolates changed from 9 (2024) to 2 (2025); possible change in panel, laboratory participation or sampling |
| WARNING | Klebsiella pneumoniae | Meropenem | Exampleland | aquatic | small_denominator_years | years with fewer than 10 isolates: 2022 (n=3), 2023 (n=6), 2024 (n=9), 2025 (n=2) |
| WARNING | Klebsiella pneumoniae | Meropenem | Exampleland | aquatic | testing_volume_jump | tested isolates changed from 3 (2022) to 6 (2023); possible change in panel, laboratory participation or sampling |
| WARNING | Klebsiella pneumoniae | Meropenem | Exampleland | aquatic | testing_volume_jump | tested isolates changed from 9 (2024) to 2 (2025); possible change in panel, laboratory participation or sampling |
| WARNING | Escherichia coli | Cefotaxime | Exampleland | human | small_denominator_years | years with fewer than 10 isolates: 2022 (n=9), 2023 (n=6), 2024 (n=6), 2025 (n=9) |
| WARNING | Escherichia coli | Cefotaxime | Exampleland | human | ast_method_changes_across_years | AST method differs between years: 2022: broth_microdilution; 2023: broth_microdilution; 2024: broth_microdilution; 2025: disk_diffusion |
| WARNING | Escherichia coli | Ciprofloxacin | Exampleland | human | small_denominator_years | years with fewer than 10 isolates: 2022 (n=9), 2023 (n=6), 2024 (n=6), 2025 (n=9) |
| WARNING | Escherichia coli | Meropenem | Exampleland | human | small_denominator_years | years with fewer than 10 isolates: 2022 (n=9), 2023 (n=6), 2024 (n=6), 2025 (n=9) |
| WARNING | Klebsiella pneumoniae | Cefotaxime | Exampleland | human | small_denominator_years | years with fewer than 10 isolates: 2024 (n=8) |
| WARNING | Klebsiella pneumoniae | Cefotaxime | Exampleland | human | testing_volume_jump | tested isolates changed from 8 (2024) to 20 (2025); possible change in panel, laboratory participation or sampling |
| WARNING | Klebsiella pneumoniae | Ciprofloxacin | Exampleland | human | small_denominator_years | years with fewer than 10 isolates: 2024 (n=8) |
| WARNING | Klebsiella pneumoniae | Ciprofloxacin | Exampleland | human | testing_volume_jump | tested isolates changed from 8 (2024) to 20 (2025); possible change in panel, laboratory participation or sampling |
| WARNING | Klebsiella pneumoniae | Meropenem | Exampleland | human | small_denominator_years | years with fewer than 10 isolates: 2024 (n=8) |
| WARNING | Klebsiella pneumoniae | Meropenem | Exampleland | human | testing_volume_jump | tested isolates changed from 8 (2024) to 20 (2025); possible change in panel, laboratory participation or sampling |

## AST category re-check against supplied breakpoints

Categories re-derived from MIC/zone and the **supplied** breakpoint rows (never built in): concordant **360**.

## Genotype–phenotype analysis

**Go / no-go decision: `GO`** — 120 isolate(s) with a single screening record and a valid AST result.

| Criterion | Met | Detail |
|---|:---:|---|
| WGS screening and AST share the same isolate (exact sample_id) | ✅ | 120 isolate(s) have a single screening record and at least one valid S/I/R result |
| AST method recorded for every paired AST record | ✅ | method must be populated |
| Breakpoint standard and version recorded for every paired AST record | ✅ | breakpoint_standard and breakpoint_version must be populated |
| Screening tool, database and database version recorded for every paired isolate | ✅ | screening table needs method, database and database_version |
| No organism conflicts between manifest and AST among paired isolates | ✅ | 0 conflicting isolate/organism label(s) |
| A gene-to-drug map with a source citation for every row was supplied | ✅ | supply --gene-drug-map with mapping_source for every row |
| Sample-selection process documented by the analyst (cannot be verified by software) | ✅ | set study.sample_selection_documented: true and describe it in study.selection_description |

Descriptive concordance of a mapped-determinant screen with recorded phenotype among **screened** isolates. Not diagnostic accuracy or causal attribution. Strata below the minimum isolate count show counts only.

| Organism | Drug | Method | Std | Ver | n | G+/R | G+/notR | G−/R | G−/notR | Sens | Spec | PPV | NPV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Escherichia coli | Cefotaxime | broth_microdilution | SYNTHETIC | 0.0 | 28 | 7 | 1 | 0 | 20 | 100.0% | 95.2% | 87.5% | 100.0% |
| Escherichia coli | Cefotaxime | disk_diffusion | SYNTHETIC | 0.0 | 12 | 4 | 0 | 0 | 8 | 100.0% | 100.0% | 100.0% | 100.0% |
| Escherichia coli | Ciprofloxacin | broth_microdilution | SYNTHETIC | 0.0 | 40 | 10 | 1 | 2 | 27 | 83.3% | 96.4% | 90.9% | 93.1% |
| Escherichia coli | Meropenem | broth_microdilution | SYNTHETIC | 0.0 | 40 | 9 | 1 | 3 | 27 | 75.0% | 96.4% | 90.0% | 90.0% |
| Klebsiella pneumoniae | Cefotaxime | broth_microdilution | SYNTHETIC | 0.0 | 80 | 15 | 2 | 2 | 61 | 88.2% | 96.8% | 88.2% | 96.8% |
| Klebsiella pneumoniae | Ciprofloxacin | broth_microdilution | SYNTHETIC | 0.0 | 80 | 27 | 1 | 3 | 49 | 90.0% | 98.0% | 96.4% | 94.2% |
| Klebsiella pneumoniae | Meropenem | broth_microdilution | SYNTHETIC | 0.0 | 80 | 24 | 2 | 3 | 51 | 88.9% | 96.2% | 92.3% | 94.4% |

21 discordant isolate–drug record(s) are listed in `genotype_phenotype_discordant.csv` for review (possible causes: unmapped mechanism, expression/regulation, breakpoint or method issue, mislabelled sample, gene not functional).

## Model predictions versus observed AST

Predictions stay a separate stream. AMR-Surveil cannot verify that these isolates were absent from model training, so accuracy-type values appear only when the analyst attests independence **and** the rows are from a held-out split.

| Model | Drug | Split | Pairs | P-R/O-R | P-R/O-notR | P-S/O-R | P-S/O-notR | Agreement | Status |
|---|---|---|---|---|---|---|---|---|---|
| EpiResNet | Cefotaxime | test | 20 | 2 | 3 | 2 | 13 | — | counts_only: independence from training data not attested or split not held-out |
| EpiResNet | Ciprofloxacin | test | 20 | 4 | 1 | 4 | 11 | — | counts_only: independence from training data not attested or split not held-out |
| EpiResNet | Meropenem | test | 20 | 5 | 1 | 1 | 13 | — | counts_only: independence from training data not attested or split not held-out |

## Missingness

Every AST record has a valid S/I/R result.

## Separate surveillance streams

- **Aggregate resistance:** 3 row(s); percentages recomputed from source counts; never joined to isolates.
- **Population incidence:** 1 row(s); source-defined cases per stated population; resistance proportion ≠ incidence.
- **Antimicrobial consumption:** 1 row(s); DDD per 1,000 inhabitants per day; not patient-level exposure.

## Interpretation checklist

- Confirm licence, accession provenance, retrieval date and the sample-selection process.
- Verify metadata and AST belong to the same isolate; record method and breakpoint standard/version.
- Review laboratory participation, specimen mix, missingness, deduplication, panel changes and breakpoint revisions before comparing years.
- Do not compare incompatible human, animal, food, aquatic or environmental denominators as one sampling frame.
- A negative gene call is not susceptibility; genomic similarity is not transmission.
- Review all QA findings and `run_manifest.json` (it hashes every input and output) before communicating results.

---
Generated by AMR-Surveil. This report describes the supplied data and configuration; it does not certify data quality, clinical validity, laboratory performance or operational readiness.
