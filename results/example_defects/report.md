# AMR-Surveil analysis report

**Software version:** `1.1.1`  
**Run timestamp (UTC):** `2026-11-14T12:19:59.957829+00:00`  
**Reference date for date checks:** `2026-11-07`  
**Deduplication rule:** `none`  
**Quality gate:** `exploratory (errors do not stop the run)`

> **Scope:** research and portfolio analysis of supplied data. Not a clinical decision-support system, a diagnostic report, a national surveillance submission, or evidence of operational validation.

> ⚠️ **Exploratory run with 4 unresolved ERROR finding(s).** Do not communicate these results until the errors are resolved.

## Study definition

- **Name:** AMR-Surveil synthetic demonstration
- **Question:** Does the pipeline recover the planted gene-phenotype association and flag the planted method change?
- **Inclusion criteria:** All synthetic isolates.
- **Exclusion criteria:** None.
- **Sample-selection description:** Synthetic data generated with a fixed seed; no real isolates or sampling frame exist.

## Dataset inventory

| Stream | Rows | Interpretation |
|---|---:|---|
| Isolate/sample manifest | 121 | Isolate-level metadata |
| Screened samples (genotyping declared) | 120 | Only these can have a *negative* genotype |
| AMR determinant calls (supplied / used) | 169 / 164 | Positive calls; evidence tiers excluded: `heuristic_motif, homology_search` |
| Phenotypic AST records | 361 | Observed categories supplied by the user |
| Aggregate surveillance rows | 4 | Source counts; not isolate-linked |
| Population incidence rows | 1 | Cases per stated denominator; not isolate-linked |
| Antimicrobial consumption rows | 1 | DDD-based indicator; not patient-level use |
| Model predictions | 60 | Research outputs; not observed AST |
| Plasmid/replicon rows | 7 | Upstream evidence; not proof of transfer |
| Isolates with AMR calls and AST (exact ID) | 87 | Co-occurrence inventory only |
| Records excluded with a recorded reason | 11 | See `exclusion_log.csv` |

## Data-quality review

Errors **4** · Warnings **8** · Info **0**

| Severity | Table | Code | Rows | Finding |
|---|---|---|---|---|
| ERROR | samples | duplicate_sample_id | 2 | sample_id must be unique in samples.csv; joins would be ambiguous. |
| ERROR | samples | future_collection_date | 1 | Collection date is after the reference date 2026-11-07. |
| ERROR | ast | invalid_ast_category | 1 | AST result must be S, I or R. EUCAST I means susceptible, increased exposure; other standards may define categories differently. |
| ERROR | aggregate_surveillance | duplicate_aggregate_key | 2 | The same aggregate-surveillance key appears more than once. Rows are never summed automatically because duplicates would double-count; remove the duplicate or add a field that distinguishes the rows. |
| WARNING | samples | inconsistent_organism_labels | 1 | Labels in 'organism' differ only by case/whitespace; they are merged onto the most frequent spelling for analysis. Synonyms (e.g. abbreviations) are not merged: curate a controlled vocabulary upstream. |
| WARNING | samples | non_iso_collection_date | 1 | Collection dates must be ISO 8601 (YYYY, YYYY-MM or YYYY-MM-DD). Other forms such as 05/06/2024 are ambiguous (day/month order) and are not parsed; the affected year is unknown. |
| WARNING | ast | unmatched_sample_id | 1 | Rows reference sample IDs absent from samples.csv; they are excluded from isolate-linked outputs. |
| WARNING | amr_calls | low_identity_or_coverage_review | 1 | Calls below identity 90% or coverage 60% are flagged for expert review. These are review flags, not call filters; no call is dropped for this reason. |
| WARNING | amr_calls | low_confidence_evidence_tier | 5 | Calls from heuristic-motif or generic homology-search methods are not validated resistance-gene detections. They are excluded from linkage and concordance by the configured evidence-tier policy. |
| WARNING | ast | missing_result | 1 | AST rows without result cannot contribute to categorical indicators. |
| WARNING | ast | invalid_mic | 1 | MIC must be a number, optionally with a censoring operator (<=, <, >, >=, =); decimal commas are not accepted. Original rows are retained for review. |
| WARNING | ast | interpretation_before_collection | 3 | An AST result is dated before its isolate's collection date. |

### Label variants merged for analysis

Case/whitespace variants were folded onto the most frequent spelling. Synonyms (e.g. *E. coli* vs *Escherichia coli*) are **not** merged.

| Table | Column | Canonical | Variants merged |
|---|---|---|---|
| samples | organism | Escherichia coli | Escherichia coli; escherichia coli |

### Exclusion summary

| Stage | Table | Reason | Records |
|---|---|---|---|
| evidence_tier | amr_calls | excluded_evidence_tier | 5 |
| indicators | ast | ambiguous_sample_id | 1 |
| indicators | ast | missing_or_invalid_result | 2 |
| indicators | ast | unmatched_sample_id | 1 |
| plasmids | plasmids | unmatched_or_ambiguous_sample_id | 1 |
| prediction_eval | predictions | unmatched_sample_id | 1 |

## Isolate-level AST resistance indicators

Resistance proportion is **R / (S + I + R)** per group. `I` is kept in the denominator and not in the numerator (EUCAST: *susceptible, increased exposure*). Wilson score intervals. The small-n flag is a display warning, not a protocol rule.

| Year | Country | Sector | Organism | Antimicrobial | Breakpoint | R / tested | Resistance (95% CI) | Caution |
|---|---|---|---|---|---|---|---|---|
| 2025 | Exampleland | human | Klebsiella pneumoniae | Meropenem | SYNTHETIC 0.0 | 7/19 | 36.8% (19.1–59.0) | — |
| 2025 | Exampleland | human | Klebsiella pneumoniae | Ciprofloxacin | SYNTHETIC 0.0 | 4/19 | 21.1% (8.5–43.3) | — |
| 2025 | Exampleland | human | Klebsiella pneumoniae | Cefotaxime | SYNTHETIC 0.0 | 4/19 | 21.1% (8.5–43.3) | — |
| 2022 | Exampleland | human | Klebsiella pneumoniae | Meropenem | SYNTHETIC 0.0 | 7/17 | 41.2% (21.6–64.0) | — |
| 2022 | Exampleland | human | Klebsiella pneumoniae | Ciprofloxacin | SYNTHETIC 0.0 | 9/17 | 52.9% (31.0–73.8) | — |
| 2022 | Exampleland | human | Klebsiella pneumoniae | Cefotaxime | SYNTHETIC 0.0 | 4/17 | 23.5% (9.6–47.3) | — |
| 2023 | Exampleland | human | Klebsiella pneumoniae | Ciprofloxacin | SYNTHETIC 0.0 | 4/14 | 28.6% (11.7–54.6) | — |
| 2023 | Exampleland | human | Klebsiella pneumoniae | Cefotaxime | SYNTHETIC 0.0 | 3/14 | 21.4% (7.6–47.6) | — |
| 2023 | Exampleland | human | Klebsiella pneumoniae | Meropenem | SYNTHETIC 0.0 | 6/14 | 42.9% (21.4–67.4) | — |
| 2024 | Exampleland | aquatic | Klebsiella pneumoniae | Meropenem | SYNTHETIC 0.0 | 2/9 | 22.2% (6.3–54.7) | Small n |
| 2022 | Exampleland | human | Escherichia coli | Cefotaxime | SYNTHETIC 0.0 | 3/9 | 33.3% (12.1–64.6) | Small n |
| 2024 | Exampleland | aquatic | Klebsiella pneumoniae | Ciprofloxacin | SYNTHETIC 0.0 | 5/9 | 55.6% (26.7–81.1) | Small n |
| 2025 | Exampleland | human | Escherichia coli | Cefotaxime | SYNTHETIC 0.0 | 3/9 | 33.3% (12.1–64.6) | Small n |
| 2025 | Exampleland | human | Escherichia coli | Meropenem | SYNTHETIC 0.0 | 1/9 | 11.1% (2.0–43.5) | Small n |
| 2024 | Exampleland | aquatic | Klebsiella pneumoniae | Cefotaxime | SYNTHETIC 0.0 | 4/9 | 44.4% (18.9–73.3) | Small n |
| 2025 | Exampleland | human | Escherichia coli | Ciprofloxacin | SYNTHETIC 0.0 | 3/9 | 33.3% (12.1–64.6) | Small n |
| 2022 | Exampleland | human | Escherichia coli | Ciprofloxacin | SYNTHETIC 0.0 | 2/9 | 22.2% (6.3–54.7) | Small n |
| 2022 | Exampleland | human | Escherichia coli | Meropenem | SYNTHETIC 0.0 | 5/9 | 55.6% (26.7–81.1) | Small n |
| 2024 | Exampleland | human | Klebsiella pneumoniae | Meropenem | SYNTHETIC 0.0 | 3/8 | 37.5% (13.7–69.4) | Small n |
| 2024 | Exampleland | human | Klebsiella pneumoniae | Ciprofloxacin | SYNTHETIC 0.0 | 1/8 | 12.5% (2.2–47.1) | Small n |
| 2024 | Exampleland | human | Klebsiella pneumoniae | Cefotaxime | SYNTHETIC 0.0 | 1/8 | 12.5% (2.2–47.1) | Small n |
| 2023 | Exampleland | human | Escherichia coli | Cefotaxime | SYNTHETIC 0.0 | 2/6 | 33.3% (9.7–70.0) | Small n |
| 2024 | Exampleland | human | Escherichia coli | Cefotaxime | SYNTHETIC 0.0 | 1/6 | 16.7% (3.0–56.4) | Small n |
| 2023 | Exampleland | aquatic | Klebsiella pneumoniae | Meropenem | SYNTHETIC 0.0 | 2/6 | 33.3% (9.7–70.0) | Small n |
| 2023 | Exampleland | aquatic | Klebsiella pneumoniae | Cefotaxime | SYNTHETIC 0.0 | 0/6 | 0.0% (0.0–39.0) | Small n |

*Showing 25 of 54 rows; see the CSV output for all rows.*

## Comparability flags

Flags mean *do not read a trend into these numbers until the cause is understood*. Nothing is adjusted or smoothed.

| Severity | Organism | Antimicrobial | Country | Sector | Flag | Detail |
|---|---|---|---|---|---|---|
| WARNING | Escherichia coli | Cefotaxime | Exampleland | aquatic | small_denominator_years | years with fewer than 10 isolates: 2022 (n=2), 2023 (n=2), 2024 (n=2), 2025 (n=3) |
| WARNING | Escherichia coli | Cefotaxime | Exampleland | aquatic | ast_method_changes_across_years | AST method differs between years: 2022: broth_microdilution; 2023: broth_microdilution; 2024: broth_microdilution; 2025: disk_diffusion |
| WARNING | Escherichia coli | Ciprofloxacin | Exampleland | aquatic | small_denominator_years | years with fewer than 10 isolates: 2022 (n=2), 2023 (n=2), 2024 (n=2), 2025 (n=3) |
| WARNING | Escherichia coli | Meropenem | Exampleland | aquatic | small_denominator_years | years with fewer than 10 isolates: 2022 (n=2), 2023 (n=2), 2024 (n=2), 2025 (n=3) |
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
| WARNING | Klebsiella pneumoniae | Cefotaxime | Exampleland | human | unknown_year_records | 1 record(s) lack a valid collection year and cannot be placed in a trend |
| WARNING | Klebsiella pneumoniae | Cefotaxime | Exampleland | human | small_denominator_years | years with fewer than 10 isolates: 2024 (n=8), 2099 (n=1) |
| WARNING | Klebsiella pneumoniae | Cefotaxime | Exampleland | human | missing_years | no isolates in year(s): 2026, 2027, 2028, 2029, 2030, 2031, 2032, 2033, 2034, 2035, 2036, 2037, 2038, 2039, 2040, 2041, 2042, 2043, 2044, 2045, 2046, 2047, 2048, 2049, 2050, 2051, 2052, 2053, 2054, 2055, 2056, 2057, 2058, 2059, 2060, 2061, 2062, 2063, 2064, 2065, 2066, 2067, 2068, 2069, 2070, 2071, 2072, 2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083, 2084, 2085, 2086, 2087, 2088, 2089, 2090, 2091, 2092, 2093, 2094, 2095, 2096, 2097, 2098 |
| WARNING | Klebsiella pneumoniae | Cefotaxime | Exampleland | human | testing_volume_jump | tested isolates changed from 8 (2024) to 19 (2025); possible change in panel, laboratory participation or sampling |
| WARNING | Klebsiella pneumoniae | Cefotaxime | Exampleland | human | testing_laboratory_set_changes | testing laboratory set differs between years: 2022: Synthetic Lab A, Synthetic Lab B; 2023: Synthetic Lab A, Synthetic Lab B; 2024: Synthetic Lab A, Synthetic Lab B; 2025: Synthetic Lab A, Synthetic Lab B; 2099: Synthetic Lab A |
| WARNING | Klebsiella pneumoniae | Cefotaxime | Exampleland | human | specimen_mix_changes | specimen types differs between years: 2022: blood, urine; 2023: blood, urine; 2024: blood, urine; 2025: blood, urine; 2099: urine |
| WARNING | Klebsiella pneumoniae | Ciprofloxacin | Exampleland | human | unknown_year_records | 1 record(s) lack a valid collection year and cannot be placed in a trend |
| WARNING | Klebsiella pneumoniae | Ciprofloxacin | Exampleland | human | small_denominator_years | years with fewer than 10 isolates: 2024 (n=8), 2099 (n=1) |
| WARNING | Klebsiella pneumoniae | Ciprofloxacin | Exampleland | human | missing_years | no isolates in year(s): 2026, 2027, 2028, 2029, 2030, 2031, 2032, 2033, 2034, 2035, 2036, 2037, 2038, 2039, 2040, 2041, 2042, 2043, 2044, 2045, 2046, 2047, 2048, 2049, 2050, 2051, 2052, 2053, 2054, 2055, 2056, 2057, 2058, 2059, 2060, 2061, 2062, 2063, 2064, 2065, 2066, 2067, 2068, 2069, 2070, 2071, 2072, 2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083, 2084, 2085, 2086, 2087, 2088, 2089, 2090, 2091, 2092, 2093, 2094, 2095, 2096, 2097, 2098 |
| WARNING | Klebsiella pneumoniae | Ciprofloxacin | Exampleland | human | testing_volume_jump | tested isolates changed from 8 (2024) to 19 (2025); possible change in panel, laboratory participation or sampling |
| WARNING | Klebsiella pneumoniae | Ciprofloxacin | Exampleland | human | testing_laboratory_set_changes | testing laboratory set differs between years: 2022: Synthetic Lab A, Synthetic Lab B; 2023: Synthetic Lab A, Synthetic Lab B; 2024: Synthetic Lab A, Synthetic Lab B; 2025: Synthetic Lab A, Synthetic Lab B; 2099: Synthetic Lab A |
| WARNING | Klebsiella pneumoniae | Ciprofloxacin | Exampleland | human | specimen_mix_changes | specimen types differs between years: 2022: blood, urine; 2023: blood, urine; 2024: blood, urine; 2025: blood, urine; 2099: urine |
| WARNING | Klebsiella pneumoniae | Meropenem | Exampleland | human | unknown_year_records | 1 record(s) lack a valid collection year and cannot be placed in a trend |

*Showing 30 of 35 rows; see the CSV output for all rows.*

## AST category re-check against supplied breakpoints

Categories re-derived from MIC/zone and the **supplied** breakpoint rows (never built in): concordant **358**, discordant **2**, no_measurement **1**.

Discordant records:

| Sample | Antimicrobial | Recorded | Derived | MIC | Zone | Standard | Version |
|---|---|---|---|---|---|---|---|
| SYN-00001 | Cefotaxime | X | S | 0.5 | — | SYNTHETIC | 0.0 |
| SYN-00001 | Ciprofloxacin | — | S | 0.25 | — | SYNTHETIC | 0.0 |

## Genotype–phenotype analysis

**Go / no-go decision: `GO`** — 119 isolate(s) with a single screening record and a valid AST result.

| Criterion | Met | Detail |
|---|:---:|---|
| WGS screening and AST share the same isolate (exact sample_id) | ✅ | 119 isolate(s) have a single screening record and at least one valid S/I/R result |
| AST method recorded for every paired AST record | ✅ | method must be populated |
| Breakpoint standard and version recorded for every paired AST record | ✅ | breakpoint_standard and breakpoint_version must be populated |
| Screening tool, database and database version recorded for every paired isolate | ✅ | screening table needs method, database and database_version |
| No organism conflicts between manifest and AST among paired isolates | ✅ | 0 conflicting isolate/organism label(s) |
| A gene-to-drug map with a source citation for every row was supplied | ✅ | supply --gene-drug-map with mapping_source for every row |
| Sample-selection process documented by the analyst (cannot be verified by software) | ✅ | set study.sample_selection_documented: true and describe it in study.selection_description |

Descriptive concordance of a mapped-determinant screen with recorded phenotype among **screened** isolates. Not diagnostic accuracy or causal attribution. Strata below the minimum isolate count show counts only.

| Organism | Drug | Method | Std | Ver | n | G+/R | G+/notR | G−/R | G−/notR | Sens | Spec | PPV | NPV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Escherichia coli | Cefotaxime | broth_microdilution | SYNTHETIC | 0.0 | 27 | 7 | 1 | 0 | 19 | 100.0% | 95.0% | 87.5% | 100.0% |
| Escherichia coli | Cefotaxime | disk_diffusion | SYNTHETIC | 0.0 | 12 | 4 | 0 | 0 | 8 | 100.0% | 100.0% | 100.0% | 100.0% |
| Escherichia coli | Ciprofloxacin | broth_microdilution | SYNTHETIC | 0.0 | 39 | 10 | 1 | 2 | 26 | 83.3% | 96.3% | 90.9% | 92.9% |
| Escherichia coli | Meropenem | broth_microdilution | SYNTHETIC | 0.0 | 39 | 9 | 1 | 3 | 26 | 75.0% | 96.3% | 90.0% | 89.7% |
| Klebsiella pneumoniae | Cefotaxime | broth_microdilution | SYNTHETIC | 0.0 | 80 | 15 | 2 | 2 | 61 | 88.2% | 96.8% | 88.2% | 96.8% |
| Klebsiella pneumoniae | Ciprofloxacin | broth_microdilution | SYNTHETIC | 0.0 | 80 | 27 | 1 | 3 | 49 | 90.0% | 98.0% | 96.4% | 94.2% |
| Klebsiella pneumoniae | Meropenem | broth_microdilution | SYNTHETIC | 0.0 | 80 | 24 | 2 | 3 | 51 | 88.9% | 96.2% | 92.3% | 94.4% |

21 discordant isolate–drug record(s) are listed in `genotype_phenotype_discordant.csv` for review (possible causes: unmapped mechanism, expression/regulation, breakpoint or method issue, mislabelled sample, gene not functional).

## Model predictions versus observed AST

Predictions stay a separate stream. AMR-Surveil cannot verify that these isolates were absent from model training, so accuracy-type values appear only when the analyst attests independence **and** the rows are from a held-out split.

| Model | Drug | Split | Pairs | P-R/O-R | P-R/O-notR | P-S/O-R | P-S/O-notR | Agreement | Status |
|---|---|---|---|---|---|---|---|---|---|
| EpiResNet | Cefotaxime | test | 19 | 2 | 3 | 2 | 12 | — | counts_only: independence from training data not attested or split not held-out |
| EpiResNet | Ciprofloxacin | test | 20 | 4 | 1 | 4 | 11 | — | counts_only: independence from training data not attested or split not held-out |
| EpiResNet | Meropenem | test | 20 | 5 | 1 | 1 | 13 | — | counts_only: independence from training data not attested or split not held-out |

## Missingness

| Country | Organism | Antimicrobial | Year | Lab | Records | Missing/invalid | % |
|---|---|---|---|---|---|---|---|
| Unknown | Escherichia coli | Cefotaxime | — | Synthetic Lab A | 1 | 1 | 100.0 |
| Unknown | Escherichia coli | Ciprofloxacin | — | Synthetic Lab A | 1 | 1 | 100.0 |

## Separate surveillance streams

- **Aggregate resistance:** 4 row(s); percentages recomputed from source counts; never joined to isolates.
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
