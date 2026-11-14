# Changelog

All notable changes to this project are documented here. Format: [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

Repository hygiene after the first public upload (no code behaviour change).
- Hidden files that a browser drag-and-drop upload silently skips were added: `.github/workflows/ci.yml`, `.github/dependabot.yml`,
  `.github/ISSUE_TEMPLATE/bug_report.md`, `.streamlit/config.toml` (the Dockerfile copies it), `.gitignore`, `.dockerignore`, `.gitattributes`.
- Code formatted with `ruff format` so the CI lint job passes.
- CI: legs not yet verified by the author (Python 3.10-3.12 and pandas 2.2) are marked experimental and do not turn the run red.
- `CITATION.cff`: removed an invalid top-level `affiliation` key, added the ORCID iD and DOI, corrected the release date.
- README: the FishAMR-Link and EpiResNet links point to their own repositories.
- CI lint job pins the ruff series used for formatting; Dependabot is quieter (no Python-image or minimum-version bumps).

## [1.1.1] - 2026-11-14

Theme: **final sweep after a second real local run**. No interface changes.

### Fixed
- **Web app: the AST re-check tab crashed with `ModuleNotFoundError: No module named 'rpds.rpds'`** on a Python 3.14 / Windows
  environment. `st.bar_chart` imports altair → jsonschema → referencing → rpds-py, and the rpds-py binary wheel was broken. The tab now
  uses a Plotly bar chart, so the app no longer needs the altair stack. A regression test forbids Streamlit's altair-backed charts,
  `amr-surveil doctor` checks `rpds`, and the user guide documents the fix (`pip uninstall -y rpds-py && pip install -U rpds-py`).
- **A blank-string `mic_unit` was treated as an unsupported unit** (`unit_not_mg_per_l`) although the documentation says blank means mg/L.
- **Wilson interval boundaries** carried floating-point residue (`0/10` gave a lower bound of 2.8e-15 instead of exactly 0).
- **CSV formula injection**: output CSVs echo user-supplied text (sample IDs, organisms, notes). Cells beginning with `=`, `+`, `@`,
  tab/CR or `-` followed by a non-digit are now prefixed with `'` so spreadsheets do not execute them; numeric columns are untouched.

### Documentation
- `DATA_DICTIONARY.md`: how to enter CLSI breakpoints (`R ≥ y` / `R ≤ y` → `mic_r_gt` / `zone_r_lt`) in the EUCAST-style columns.
- `REFERENCES.md`: added Wilson (1927), CLSI M39 (5th ed.), AMRFinderPlus (Feldgarden 2021) and BV-BRC (Olson 2023); DOIs checked.
- `INTEGRATION_GUIDE.md`: why AMR-UNBIND (steered-MD descriptors for one protein–ligand pair) is deliberately not ingested.
- 250 further seeded damaged-input runs (new seeds): zero unexpected crashes.

## [1.1.0] - 2026-11-12

Theme: **first-run fixes from a real local session, and public-data access**. Everything is additive; the 1.0.0 contract is unchanged
(`antibiogram.csv` and the `doctor` / `public-data` commands are new, minor-version additions).

### Fixed (reported from a real Windows / VS Code / conda run of 1.0.0)
- **The Comparability and Data-quality tabs crashed** with `SyntaxError: '(' was never closed`: Streamlit "magic" re-parses a bare
  conditional-expression statement (`st.dataframe(x) if cond else st.success(y)`). Replaced by explicit `if/else` (`show_table`).
  A static test now scans `app.py` for such statements so the stub-based UI tests' blind spot is closed.
- `use_container_width` deprecation warnings: a compatibility helper uses `width="stretch"` and falls back on older Streamlit.
- Year axis showed fractional ticks ("2,022.5"): the axis is now categorical.
- The local config bound the server to `0.0.0.0` (reachable from the network); it now binds `127.0.0.1`. The Docker image still binds all interfaces explicitly.
- A broken `pyarrow` wheel (Python 3.14 / conda+pip mix) produced an opaque stack trace: `amr-surveil doctor` now diagnoses it and prints the fix.

### Added
- **Public data hub** (`amr_surveil.datasources`, app tab and `amr-surveil public-data ...`): live BV-BRC query (bounded), BV-BRC / NCBI Pathogen Detection / ECDC Atlas file import, converted to canonical inputs with a provenance JSON (query URL, retrieval time, response hash) and the source's terms written to `data_license`. Fetching is HTTPS-only, host-allow-listed, size-bounded, cached with hash re-check. Nothing is bundled.
- **Map tab**: per-country resistance for one organism / drug / sector / test-method context; counts summed (never percentages averaged), countries below a minimum isolate count suppressed, unrecognised country names reported.
- **Cumulative antibiogram** (CLSI M39 style, % susceptible, minimum 30 isolates) as a tab and `antibiogram.csv` output.
- **Help tab** and [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md); [`docs/RELATED_WORK.md`](docs/RELATED_WORK.md).
- `amr-surveil doctor` (environment check plus a 24-isolate smoke run).

### Fixed in the new code before release (independent adversarial review of the public-data modules, 11 defects)
- Antibiogram tab crashed on undated data (typical BV-BRC export); host sectors were mislabelled ("non-human primate" → human, "Humans" → animal, placeholders → animal); NCBI rows with a duplicate identifier leaked AST into the kept isolate and parsing was quadratic; country spellings of one country were suppressed separately on the map; BV-BRC country/year filters ignored aliases or dropped records silently; ECDC thousands separators, decimals and conflicting duplicates were lost or truncated silently; truncated downloads escaped as raw exceptions; a damaged cache entry was never repaired and the cache directory was shared and predictable; ISO alias gaps ("Sudan (South)" resolved to Sudan). Each has a regression test. The antibiogram now keeps sectors in separate rows.

### Changed
- Removed `docs/HANDBOOK_TRACEABILITY.md` and "handbook §" citations: they referred to a private learning plan readers cannot see. Rules now cite what they implement.

### Not verified
- The public-data connectors were tested with mocked transport only; the sandbox that built them had no access to the live servers. Field and column names follow each provider's published documentation (see `docs/LIMITATIONS.md`).

## [1.0.0] - 2026-11-11

Theme: **stable contract and hardening**. 1.0.0 means the CLI, the configuration keys, the input/output file names and columns,
the exclusion-reason vocabulary and the manifest format are now a public, versioned contract ([`docs/STABILITY.md`](docs/STABILITY.md)).
It does **not** mean validated on real surveillance data (see [`docs/VALIDATION.md`](docs/VALIDATION.md)).

### Added
- `amr-surveil verify RESULTS_DIR [--check-inputs]`: re-hash every recorded output (and optionally input) against `run_manifest.json`; exit 1 on any mismatch.
- `amr-surveil templates --out DIR`: header-only input templates for every table plus a default `config.yaml` (never overwrites without `--force`).
- `manifest_schema_version` in `run_manifest.json`; `n_excluded_uncertain_genotype` in the concordance table; AST re-check status `unit_not_mg_per_l`.
- Validation: `non_finite_*`, `cases_exceed_population`, `implausible_year`, `implausible_period_days`, partial-date `interpretation_date`.
- Tests: independent stdlib recomputation of indicators/concordance (`test_independent_oracle.py`), seeded fuzzing of damaged inputs (`test_robustness.py`), documented CLI commands must parse (`test_docs_sync.py`), regression tests for every defect below.

### Fixed (found by fuzzing and by two independent adversarial review rounds)
- A record without `sample_id` crashed plasmid filtering under pandas 3; such records are now logged (`missing_sample_id`) and removed.
- Deduplication ignored whether an isolate was tested: an untested first isolate removed a tested repeat. It now ranks only isolates with an eligible AST record; a year-/month-only date makes "first" uncertain and sets `date_tie`; the rule now also applies to the genotype–phenotype analysis.
- One species split into two indicator groups when AST and manifest spelled the organism differently; labels are now harmonised jointly across tables.
- A MIC in a non-mg/L unit was compared with mg/L breakpoints; a record with missing standard/version matched a breakpoint row with missing standard/version; a blank AST organism was not filled from the manifest before the re-check.
- Isolates whose only mapped gene hit was weak-tier or AMRFinderPlus `PARTIAL`/`INTERNAL_STOP` were counted as genotype-**negative** (false negatives, deflated sensitivity); they are now excluded and counted separately.
- Predictions: conflicting observed AST, unrecognised `predicted_class` and unmatched IDs were dropped silently (now logged); train and test predictions of the same isolate/drug were treated as duplicates.
- A blank antimicrobial was pooled as "Unspecified"; now excluded and logged.
- Evidence tiers were matched case-sensitively in the pipeline but case-insensitively in the audit (upper-case tiers leaked through the filter).
- Collection/interpretation dates outside 1677–2262 crashed `audit`/`qa` with `OutOfBoundsDatetime`; they are now reported findings.
- `duplicate_consumption_key` ignored `period_days` (false stop on monthly rows); infinite values passed QA and produced infinite indicators; the report always said "95% CI" whatever the configured level; invalid context-stream rows vanished without an exclusion-log entry.
- FishAMR-Link feature table with a blank `feature_type` crashed; matching is now case-insensitive. EpiResNet `binary_valid` accepts `1`/`yes`. AMRFinderPlus sample ID defaults to the file stem (so `E.coli_1.tsv` is not read as `E`); MIC `1.` now parses.

### Changed
- FishAMR-Link feature-table rows whose mechanism is `PARTIAL*`/`INTERNAL_STOP*` get evidence tier `homology_search` instead of `curated_database`.

## [0.2.0] - 2026-11-09

Theme: from a descriptive integration layer to an **auditable genotype–phenotype surveillance workbench** with
screened-sample genotypes, a go/no-go gate, comparability flags and real-format adapters.

### Added
- **Genotype–phenotype concordance** (`genotype_phenotype.py`): 2×2 per organism/drug/method/breakpoint/screening stratum with
  Wilson intervals, discordant-isolate listing, `I` counted as not-R and reported separately, small strata show counts only.
- **Go/no-go gate** with seven explicit criteria; concordance is computed only on GO.
- **Screening table**: a *negative* genotype exists only for samples declared as screened with one method/database/version.
- **Evidence tiers** for AMR calls (`curated_database`, `homology_search`, `heuristic_motif`, `unspecified`); low-evidence tiers are excluded
  by default and every exclusion is logged.
- **Comparability flags**: unknown year, single year, small-n years, gaps, testing-volume jumps, and changes in
  breakpoint version, AST method, laboratory set and specimen mix across years.
- **AST category re-check** against user-supplied breakpoints with censored-MIC logic and zone diameters.
- **Data inventory**, **AST missingness** summary and a complete **exclusion log**; label harmonisation record.
- **Prediction-vs-observed** descriptive agreement that stays counts-only unless the analyst attests independence *and* the rows are from a held-out split.
- **Adapters for real exports**: AMRFinderPlus TSV (v3.x and v4.x headers), FishAMR-Link feature table and SQLite warehouse
  (read-only; seeded simulated baseline excluded; coordinates never exported), EpiResNet `test_predictions.csv` with threshold,
  temperature, split strategy, training-data hash and checkpoint SHA-256.
- Strict ISO-8601 date parsing, MIC parsing with censoring operators, deterministic label harmonisation.
- Specimen-aware deduplication rule (`first_per_patient_species_specimen_year`); every sample gets a logged decision (kept/removed/not_assignable) and date ties are flagged.
- `study:` configuration block (name, question, criteria, sample-selection attestation) echoed into the report.
- `qa` CLI command (exit status 1 on ERROR findings, for CI), `normalize-amrfinder`, extended `normalize-fishamr` / `normalize-epiresnet`.
- Run manifest now records Python/pandas versions, the effective configuration and SHA-256 of **every output**.
- New Streamlit tabs (comparability, genotype–phenotype, AST re-check, exclusions, inventory); UI logic moved to a tested module.
- JSON Schemas generated from the code (CI fails on drift); documentation-sync tests; 218 tests including a mutation-checked core.
- Benchmarks now include a ground-truth check: the planted genotype–phenotype association must be recovered at every size.

### Fixed (defects in 0.1.0 found by adversarial probing; regression tests in `tests/test_v010_regressions.py`)
- A **blank AST result** was reported as an ERROR (`invalid_ast_category`) and stopped strict runs; it is now a `missing_result` WARNING.
- **Censored MICs** (`<=0.5`, `>32`) were flagged invalid; they are now parsed and reported as INFO.
- **Label case variants** (`Cefotaxime`, `cefotaxime`, `CEFOTAXIME `) silently split strata; they are now merged and reported.
- **Duplicate aggregate rows** double-counted silently; duplicate aggregate/incidence/consumption keys are now ERRORs.
- **Future dates** and ambiguous `dd/mm/yyyy` dates were accepted without an ERROR; future dates are ERRORs, non-ISO dates WARNINGs.
- The v0.1.0 demo README documented `--out/--samples/--seed` flags that `scripts/generate_demo.py` did not accept.
- Duplicate-key checks silently skipped when an optional key column was absent.
- `prediction` contract lacked the `organism` column that linkage already used.
- Temporary `exploratory_config.yaml` files are gone; exploratory mode is a validated in-memory override recorded in the manifest.

### Changed
- Default synthetic demo is now 120 samples with planted structure (association, 2025 method change, repeat patients, low-evidence calls).
- `pandas` range widened to `>=2.2.3,<4`; verified locally on pandas 3.0.5 only (CI matrix also targets 2.2.3).
- `ruff format` applied repository-wide; CI enforces it.

### Removed
- `docs/VALIDATION_PLAN.md` (superseded by `docs/VALIDATION.md`, which records evidence as well as the plan).

### Known limits
See [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md). In particular: Docker build and the Streamlit runtime were not executed in the authoring
environment; adapter fixtures follow documented export layouts and are not captures from live tool runs.

## [0.1.0] - 2026-11-07

### Added
- Canonical CSV data contracts for isolate manifests, AMR determinant calls, plasmid evidence, observed AST, aggregate resistance, population incidence, antimicrobial consumption and model predictions.
- Surveillance-specific QA checks, strict quality gate and machine-readable findings.
- Isolate-level R/(S+I+R) indicators with Wilson confidence intervals and retained AST method/breakpoint context.
- Exact `sample_id` evidence linkage; aggregate surveillance and model predictions remain distinct streams.
- Configurable FishAMR-Link and EpiResNet export adapters.
- Optional source-cited gene–drug descriptive stratification among gene-positive isolates.
- One Health sector context, source-specific denominator safeguards, a provenance manifest with input SHA-256 digests and deterministic synthetic fixtures.
- CLI, Streamlit dashboard, Docker/Hugging Face Space configuration, CI, tests and benchmark harness.

### Limitations
- Raw FASTQ assembly/annotation and external AMR caller execution are not included in this initial release.
- The bundled fixtures are synthetic. No clinical or operational surveillance validation is claimed.
