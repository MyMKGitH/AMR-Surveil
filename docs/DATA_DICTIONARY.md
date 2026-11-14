# Data dictionary and CSV contracts

All CSVs are UTF-8 with a header row; `sample_id` is read as a string (leading zeros are preserved);
blank cells are missing values; column names are case-sensitive after whitespace is stripped. Duplicate
or blank header names are rejected. Machine-readable JSON Schemas are in [`schemas/`](../schemas) (generated from
the code by `scripts/generate_schemas.py`; CI fails if they drift) and runnable examples are in [`data/demo/`](../data/demo).
`tests/test_docs_sync.py` fails if a column in the code is missing from this file.

**Dates** must be ISO 8601 (`YYYY`, `YYYY-MM`, `YYYY-MM-DD`). Ambiguous forms such as `05/06/2024` are never guessed:
they are reported (`non_iso_collection_date`) and the record is placed in an *unknown year* group.
**MIC** values may carry one operator (`<=`, `<`, `>=`, `>`, `≤`, `≥`); censored values are never treated as exact.
Decimal commas (`1,5`) are rejected rather than reinterpreted.

## 1. `samples.csv` — required

| Field | Required | Meaning |
|---|---|---|
| `sample_id` | **Yes** | Unique isolate/sample identifier; the only key used to link isolate-level streams. |
| `organism` | Recommended | Reported taxon. Case/whitespace variants are merged for analysis and reported; synonyms are not merged. |
| `collection_date` | Recommended | ISO date of sample collection. |
| `country`, `region` | Recommended | Geography at an ethically appropriate resolution. |
| `source_sector` | Recommended | `human`, `animal`, `veterinary`, `food`, `environment`, `aquatic`, `aquaculture`, `wastewater`, `unknown` (other labels are kept and flagged). |
| `specimen` | Recommended | Blood, urine, water, fish tissue… used by the optional specimen-aware deduplication rule and comparability checks. |
| `patient_id` | Optional | Pseudonymised key; needed only for patient-based deduplication. |
| `host` | Optional | Host or source descriptor. |
| `sequencing_accession` | Optional | Public accession for the sequence data. |
| `sequencing_type` | Optional | e.g. Illumina short-read, ONT, hybrid. |
| `source_study`, `data_license` | Recommended | Study identifier and licence/terms (record licence and retrieval provenance). |

## 2. `amr_calls.csv` — optional (positive hits only)

| Field | Required | Meaning |
|---|---|---|
| `sample_id` | **Yes** | ID in the manifest. |
| `gene` | **Yes** | Determinant symbol exactly as called. |
| `organism` | Optional | Caller-assigned taxon; conflicts with the manifest exclude the row from genotype analysis. |
| `identity_pct`, `coverage_pct` | Optional | 0–100. Values below `review_identity_pct`/`review_coverage_pct` are flagged for review, never silently dropped. |
| `contig`, `start`, `end` | Optional | Coordinates. |
| `method` | Recommended | Caller and method. |
| `database`, `database_version` | Recommended | Database and release. |
| `source_system` | Recommended | Exporting tool. |
| `element_type` | Optional | AMRFinderPlus element type (`AMR`, `STRESS`, `VIRULENCE`…). The AMRFinderPlus adapter keeps `AMR` by default. |
| `determinant_class` | Optional | Drug class/subclass reported by the caller (informational; **not** used as the gene→drug map). |
| `hit_method` | Optional | Caller's own method code (AMRFinderPlus: `EXACTX`, `ALLELEX`, `BLASTX`, `PARTIALX`, `HMM`, `POINTX`, `INTERNAL_STOP`…). Partial and internal-stop hits are flagged. |
| `evidence_tier` | Recommended | `curated_database` (alignment/HMM/point-mutation search against a curated, versioned database), `homology_search` (generic similarity, e.g. remote BLAST), `heuristic_motif` (peptide/regex motif), `unspecified`. By default `heuristic_motif` and `homology_search` are **excluded** from genotype analysis (each exclusion is logged). |

A row is a **positive** call. A missing row means nothing unless the sample is declared in `screening.csv`.

## 3. `screening.csv` — optional, required for any genotype–phenotype analysis

One row per sample that was actually screened. This is what turns "no calls" into a *negative* instead of *missing data*.

| Field | Required | Meaning |
|---|---|---|
| `sample_id` | **Yes** | Screened sample. |
| `method`, `tool_version` | Recommended | Tool and version (AMRFinderPlus does not write its version into the report: supply it). |
| `database`, `database_version` | Recommended | Database snapshot. Required by the go/no-go gate. |
| `identity_threshold_pct`, `coverage_threshold_pct` | Optional | Thresholds applied by the caller. |
| `source_system` | Optional | Exporting system. |

Samples with more than one distinct method/database/version are excluded from concordance (and logged).

## 4. `ast.csv` — optional

| Field | Required | Meaning |
|---|---|---|
| `sample_id` | **Yes** | Isolate ID. |
| `organism` | Optional | AST-record taxon; conflicts with the manifest are excluded and logged. |
| `antimicrobial` | **Yes** | Drug label (case/whitespace variants merged and reported). |
| `result` | **Yes** (column) | `S`, `I` or `R`. Blank = missing (WARNING); anything else = ERROR. EUCAST `I` = *susceptible, increased exposure*. |
| `mic`, `mic_unit` | Optional | Measured MIC with optional censoring operator. |
| `zone_mm` | Optional | Disk-diffusion zone diameter (positive number). |
| `method` | Recommended | e.g. broth_microdilution, disk_diffusion. Never pooled across methods. |
| `breakpoint_standard`, `breakpoint_version` | Recommended | e.g. EUCAST, 14.0. Never pooled across versions. |
| `interpretation_date` | Optional | ISO date; earlier than collection = WARNING. |
| `testing_lab` | Optional | Used for comparability and missingness summaries. |

## 5. `breakpoints.csv` — optional

Used only to **re-derive** categories from MIC/zone and compare them with the recorded category. No breakpoints are built in.

| Field | Meaning |
|---|---|
| `organism`, `antimicrobial`, `standard`, `version` | Join key (case-insensitive) to AST rows. Ambiguous duplicates are ignored, never guessed. |
| `mic_s_le`, `mic_r_gt` | MIC ≤ `mic_s_le` → S; > `mic_r_gt` → R; between → I. |
| `zone_s_ge`, `zone_r_lt` | Zone ≥ `zone_s_ge` → S; < `zone_r_lt` → R; between → I. |
| `source_url` | Where the values were taken from. |

**Entering CLSI breakpoints.** The columns use the EUCAST-style convention (S: MIC ≤ value; R: MIC > value; zone S ≥ value; R < value). CLSI tables are written as "S ≤ x, R ≥ y" for MICs and "S ≥ x, R ≤ y" for zones, so convert before entering them: `mic_r_gt` is the highest MIC that is still *not* resistant (the doubling dilution below `y`), and `zone_r_lt` is `y + 1` mm. Where a standard defines "I" differently (EUCAST: susceptible, increased exposure; CLSI: intermediate), the recorded category is compared as the letter `I` only.

## 6. `gene_drug_map.csv` — optional, required for concordance

| Field | Required | Meaning |
|---|---|---|
| `organism`, `gene`, `antimicrobial` | **Yes** | A gene counts as positive evidence for that drug in that organism. |
| `mapping_source` | **Yes** | Citation. Rows without one make the go/no-go gate fail. No map is bundled: curating it is a scientific decision. |
| `match_type` | Optional | `exact` (default) or `prefix` (e.g. `blaCTX-M` covers `blaCTX-M-15`). |

## 7. `plasmid_calls.csv` — optional

`sample_id`, `replicon` (required); `mobility_prediction`, `method`, `database`, `database_version`, `source_system`,
`evidence_tier` (same vocabulary as above; excluded tiers apply). A replicon is not proof of transfer.

## 8. `aggregate_surveillance.csv` — optional, separate stream

`year`, `country`, `organism`, `antimicrobial`, `resistant_n`, `tested_n`, `source_name`, `metric_definition` (required);
`source_sector`, `source_url` (optional). Percentages are recomputed from counts. Duplicate keys are an ERROR
(rows are never summed).

## 9. `incidence.csv` — optional, separate stream

`year`, `country`, `organism`, `antimicrobial`, `resistant_infection_n`, `population_n`, `source_name`, `case_definition`
(required); `source_sector`, `population_multiplier` (default 100,000), `source_url`. Zero/negative/fractional values are errors.

## 10. `consumption.csv` — optional, separate stream

`year`, `country`, `care_sector`, `antimicrobial`, `ddd_total`, `population_n`, `period_days`, `source_name` (required);
`atc_code`, `aware_group`, `data_source_type` (`sales`, `distribution`, `dispensing`, `reimbursement`, `other`), `source_url`.
DDD per 1,000 inhabitants per day = `ddd_total / population_n / period_days × 1000`.

## 11. `model_predictions.csv` — optional, separate stream

| Field | Meaning |
|---|---|
| `sample_id` (required), `organism`, `antimicrobial` (required) | Exact-ID link; organism conflicts excluded. |
| `predicted_mic`, `predicted_mic_unit` | Predicted MIC if the model produces one (EpiResNet v5 training outputs do not). |
| `predicted_class` | `predicted_resistant` / `predicted_susceptible`; anything else is retained but not scored. |
| `prediction_score` | Model probability (0–1). **Not** a statistical confidence. |
| `model_name`, `model_version`, `checkpoint_id`, `training_data_version`, `split_strategy`, `source_system` | Provenance. `checkpoint_id` is the SHA-256 of the checkpoint file in the EpiResNet adapter. |
| `prediction_unit` | e.g. `protein-antibiotic pair`: the unit the model predicts. Several predictions for one isolate/drug are **not** collapsed into an invented isolate label. |
| `protein_cluster` | Cluster ID, for leakage review. |
| `decision_threshold`, `probability_temperature` | Calibration/threshold used to produce `predicted_class`. |
| `evaluation_split` | `test`/`held_out` for rows from a held-out split. Needed (with analyst attestation) before agreement metrics are shown. |

## Outputs

| File | Content |
|---|---|
| `audit_issues.csv` | Severity-sorted QA findings (`ERROR` / `WARNING` / `INFO`). |
| `exclusion_log.csv` | Every record excluded from an analysis, with stage and reason. |
| `label_harmonisation.csv` | Case/whitespace label variants that were merged. |
| `data_inventory.csv` | One row per manifest sample: provenance and which evidence streams exist. |
| `ast_missingness.csv` | AST completeness by country, organism, drug, year, laboratory. |
| `resistance_by_group.csv` | R/(S+I+R) with Wilson CI, by year/country/sector/organism/drug/method/breakpoint. |
| `antibiogram.csv` | Cumulative antibiogram over all years: % susceptible = S/(S+I+R) per organism/drug/sector/method/breakpoint (sectors never pooled; countries pooled) with Wilson CI and a `meets_min_n` flag (default 30 isolates, CLSI M39 guidance). I counts as not susceptible here, unlike the resistance proportion. |
| `comparability_flags.csv` | Threats to cross-year comparability (see METHODS). |
| `ast_category_recheck.csv` | Recorded vs re-derived categories (when breakpoints are supplied). |
| `deduplication_log.csv` | Decision (kept / removed / not_assignable) for every sample under the chosen rule. |
| `readiness_assessment.json` | Go/no-go decision with each criterion. |
| `genotype_phenotype_concordance.csv`, `genotype_phenotype_discordant.csv` | Only when the gate says GO. Concordance columns include `n_I_in_notR` (intermediate isolates counted as not-R) and `n_excluded_uncertain_genotype` (isolates left out because their only mapped hit was weak-tier or partial/internal-stop). |
| `matched_isolates.csv`, `gene_drug_descriptive_stratification.csv` | Exact-ID co-occurrence inventory; descriptive strata among gene-positive isolates. |
| `model_prediction_context.csv`, `prediction_vs_observed.csv` | Predictions with metadata; agreement counts (metrics only if attested). |
| `plasmid_calls_linked.csv` | Replicon rows linked to unique manifest IDs. |
| `aggregate_surveillance_context.csv`, `incidence_context.csv`, `consumption_context.csv` | Separate streams with computed indicators. |
| `analysis_summary.json`, `report.md` | Machine- and human-readable summaries. |
| `run_manifest.json` | Versions, effective configuration, SHA-256 of every input **and** output. |
