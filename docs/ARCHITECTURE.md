# Architecture

```text
 FishAMR-Link ──┐                                                       ┌─ resistance_by_group.csv, comparability_flags.csv
 AMRFinderPlus ─┤ adapters ──► canonical tables ──► audit ──► gate ──► ─┤ genotype_phenotype_concordance.csv (GO only)
 EpiResNet ─────┤                 (CSV)            (strict)            ├─ prediction_vs_observed.csv, *_context.csv
 your lab files ┘                                                       └─ report.md, run_manifest.json (SHA-256 of all I/O)
```

## Principles

1. **Separate evidence types.** Genomic calls, phenotypic AST, model predictions, aggregate resistance, incidence and consumption have different denominators and meanings, so they have different tables and are never joined by country/year/organism.
2. **Say what is missing.** Absence is represented explicitly (screening table, evidence tier, `NA`), never inferred.
3. **Stop or say so.** The strict gate stops on ERROR findings. The genotype–phenotype gate refuses to compute concordance if its preconditions fail. Exploratory mode continues but records that it did.
4. **Reproducible by construction.** Seeded synthetic data, validated configuration, deterministic outputs, input and output hashes.

## Modules (`amr_surveil/`)

| Module | Role |
|---|---|
| `constants.py`, `config.py` | Column contracts, vocabularies, validated YAML configuration with documented defaults (`config/default.yaml`) |
| `io.py`, `normalize.py`, `exceptions.py` | CSV reading (strict, string-typed), label/date/MIC parsing, error hierarchy (`ConfigurationError`, `InputFileError`, `QualityGateError`) |
| `adapters/` | `amrfinder.py`, `fishamr.py`, `epiresnet.py` (+ `common.py` JSON column mapping) → canonical tables |
| `validation.py` | ~30 QA rules with severity (ERROR/WARNING/INFO) and codes |
| `exclusions.py` | `ExclusionLog`: every excluded record with its reason |
| `ast_interpretation.py` | MIC parsing with censoring; category re-check against supplied breakpoints |
| `epidemiology.py` | Eligibility, deduplication, resistance proportions, Wilson CI |
| `comparability.py` | Trend-comparability flags (§4 of METHODS) |
| `genotype_phenotype.py` | Go/no-go gate; concordance; discordant listing |
| `linkage.py`, `inventory.py` | Exact-ID co-occurrence, plasmid linkage, descriptive gene–drug strata; data inventory and missingness |
| `prediction_eval.py` | Prediction-vs-observed with attestation logic |
| `surveillance.py` | Aggregate, incidence, consumption streams |
| `pipeline.py` | Orchestration; writes 22 output files and the manifest |
| `reporting.py`, `manifest.py` | Markdown report; unique-sample-ID helpers |
| `provenance.py` | `verify` (re-hash outputs/inputs against `run_manifest.json`) and `templates` (header-only inputs + default config) |
| `synthetic.py` | Deterministic synthetic data with planted truth and optional injected defects |
| `ui_support.py`, `app.py` | Streamlit UI (logic in a tested module; `app.py` is layout only) |
| `cli.py` | `run`, `audit`, `qa`, `demo`, `verify`, `templates`, `normalize-*` |

## Pipeline steps

read → fill canonical columns → audit raw tables (strict gate) → harmonise labels → evidence-tier filter (logged) → eligibility + deduplication (logged) → indicators → comparability flags → AST re-check → inventory, missingness, linkage → prediction evaluation → plasmids → readiness gate → concordance (GO only) → write outputs, summary, report, manifest.

## Deployment

`Dockerfile` (python:3.12-slim, non-root user, port 7860, healthcheck) serves the Streamlit app and also contains the CLI; `docker-compose.yml` for local use; Hugging Face Docker Space compatible (front-matter in `README.md`). Do not upload identifiable data to a public Space.
