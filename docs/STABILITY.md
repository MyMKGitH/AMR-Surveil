# Stability contract (what "1.x" promises)

Semantic versioning applies from 1.0.0. **1.0.0 is a promise about the interface, not a claim of validation on real data** — see
[`VALIDATION.md`](VALIDATION.md) for the evidence and for what has not been done.

## Stable public surface (changes only in a major version)

| Surface | Contract |
|---|---|
| CLI | Command names `run`, `audit`, `qa`, `demo`, `verify`, `templates`, `normalize-amrfinder`, `normalize-fishamr`, `normalize-epiresnet`, and (from 1.1.0) `doctor` and `public-data`; their documented flags; exit codes `0` success, `1` QA ERROR findings (`qa`) / verification mismatch (`verify`), `2` input, configuration or quality-gate error |
| Input tables | File contracts in [`DATA_DICTIONARY.md`](DATA_DICTIONARY.md): required columns and their meaning. New *optional* columns may be added in a minor version |
| Configuration | Keys in `config/default.yaml` (`analysis`, `quality`, `study`). Unknown keys are rejected. Defaults may change only in a major version |
| Outputs | The 23 output file names (`antibiogram.csv` added in 1.1.0) and their existing column names. New columns may be appended in a minor version; none are removed or renamed |
| Vocabularies | Exclusion reasons, QA finding codes, AST re-check statuses, evidence tiers, deduplication-rule names, readiness-gate criteria. New values may be added in a minor version |
| Manifest | `run_manifest.json` with `manifest_schema_version: 1` |
| Python API | `amr_surveil.pipeline.run_pipeline`, `amr_surveil.config.load_config`, and the `amr_surveil.adapters.*` normalise functions. Other modules are internal and may change |

## What can change in minor/patch releases
Bug fixes that change a *wrong* number (these are listed in the changelog with a regression test), new QA findings, new adapters,
new optional columns, performance, documentation.

## Scientific results are versioned too
Because a bug fix can change an indicator, every run records the software version, effective configuration and SHA-256 of every
input and output; `amr-surveil verify` re-checks them. Cite the version with the results.

## Not covered by the contract
The public-data connectors' remote field and column names (they follow third-party services that can change; they are tested with mocked transport only), Streamlit page layout, report prose, the synthetic data generator's exact random draws, and anything described as a "known gap" in
[`LIMITATIONS.md`](LIMITATIONS.md).

## Recommended practice before tagging a release
Push, let CI go green on the full matrix (including the Docker build and health check), then tag. Until CI has run, treat the
build as a release candidate — those steps could not be executed where the code was written. The first public tag can be `v1.1.1`:
`1.0.0` was never published.
