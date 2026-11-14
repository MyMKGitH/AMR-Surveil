# Limitations and honest boundaries

## What this software is
A research/portfolio **integration and audit layer** for AMR surveillance evidence. It is not a clinical decision-support system, a
laboratory information system, an accredited diagnostic, a national surveillance submission tool, or a validated operational pipeline.

## What was and was not verified in v1.1.0

| Claim | Status |
|---|---|
| 388 automated tests pass with Python warnings as errors (Python 3.13.18, pandas 3.0.5) | **Verified** |
| ~96% line coverage (measured on 1.0.0; **not re-measured** for 1.1.0) | **Approximate**: measured with a stdlib tracer because `coverage.py` was not installable offline; CI enforces a floor with pytest-cov |
| Pipeline recovers a planted genotype–phenotype association at 100–20,000 isolates | **Verified** (benchmark; pooled PPV 0.92, planted 0.92) — on synthetic data only |
| Defects found in v0.1.0 are fixed | **Verified** by reproduction on the original code and regression tests |
| Core logic is exercised by tests (6 injected faults in Wilson, MIC logic, concordance counts, tier filter, small-n boundary, gate) | **Verified** (mutation spot check; not a full mutation run) |
| Docker image builds and runs | **Not verified** — no Docker daemon was available. CI contains build + smoke-test steps. |
| Streamlit dashboard renders correctly | **Partly verified.** The author ran 1.0.0 in a real browser (Windows, Python 3.14) and reported two defects, both fixed in 1.1.0 and guarded by static tests. 1.1.0 itself (new Map, Antibiogram, Public data hub tabs and `plotly` choropleth) was exercised only with stubbed `streamlit`/`plotly`; no real-browser run. |
| Public-data connectors (BV-BRC, NCBI Pathogen Detection, ECDC Atlas) | **Not verified against live services.** The build environment could not reach them. Tested: URL construction and escaping, host allow-list, size/decompression limits, caching, parsing of hand-written fixtures that follow the providers' documentation, and that converted data run through the pipeline. Column or field names (e.g. BV-BRC `genome_amr`/`genome` fields, ECDC Atlas `Indicator`/`NumValue`, NCBI `AST_phenotypes` format `drug=R,...`) may differ from the live data; failures say which field is missing. |
| pandas 2.2.x and Python 3.10–3.12 | **Not run** here. The CI matrix targets them; treat them as untested until CI is green. |
| Adapters handle real FishAMR-Link / AMRFinderPlus / EpiResNet output | **Partially verified**: fixtures follow the *documented* layouts and the upstream source code, but are hand-written, not captures from live runs. AMRFinderPlus v4 header aliases are from documentation and not tested against a v4 binary. |
| Biological validity | **Not claimed.** No real isolate data, no external benchmark, no accuracy claim. |

## Scientific boundaries

- **Negative ≠ susceptible.** A "no calls" result is a negative genotype only for a declared screen, and even then says nothing about unmapped mechanisms.
- **Gene–drug maps are not provided.** They encode scientific judgement and must be cited and reviewed by a microbiologist. A bundled map would be an unvalidated claim.
- **No breakpoints are provided.** They change between versions; supply, cite and version them.
- **Concordance inherits selection bias.** Isolates sent for sequencing are rarely representative of a surveillance population.
- **Predictions cannot be validated here.** Independence from training data cannot be checked by this tool; EpiResNet rows are protein-level, not isolate-level.
- **Genomic similarity ≠ transmission.** No typing/phylogeny is performed.
- **Streams are never joined ecologically.** Aggregate, incidence and consumption data are not linked to isolates by country/year/organism.
- **Representativeness is a design question.** Sampling frames, laboratory participation and case definitions decide what a proportion means; the flags only surface some threats.
- **Small samples:** intervals are shown, but nothing is adjusted for clustering or repeated sampling.

## Public-data boundaries

- Repository genomes are a **convenience sample** (outbreak- and resistance-enriched). Proportions are not prevalence; do not rank countries.
- NCBI/BV-BRC tables carry no patient ID: deduplication is impossible. NCBI tables usually lack AST method and breakpoint standard (`unspecified`), so comparability is unknown.
- `Nonsusceptible` phenotypes are dropped (not recoded) because they merge I and R. `Susceptible-dose dependent` is mapped to `I`.
- ECDC rows with percentages only are skipped; counts are never reconstructed from rounded percentages.
- Terms of each source were read on 2026-11-12 and can change; the app shows them next to each dataset. NCBI states it cannot grant permission on submitters' behalf, so nothing is bundled.
- Country names are matched to ISO 3166-1 alpha-3 from a built-in list; unrecognised names are reported, not guessed.

## Known gaps (reported by independent review, deliberately not "fixed" silently)

- **CLSI `SDD` / `NS` categories** are not part of the S/I/R contract and are an ERROR in strict mode; map them upstream (e.g. SDD → I) and record the rule.
- A censored MIC that is unresolved does not fall back to a disk zone, because the two come from different methods.
- `patient_id` is folded to lower case for deduplication; if your IDs are case-sensitive and differ only by case, deduplication will merge those patients.
- `matched_isolates.csv` (descriptive exact-ID co-occurrence) is built from all AST records, not the de-duplicated set.
- AMRFinderPlus rows with a blank gene symbol fall back to the free-text sequence name, which will never match a gene map.
- A TSV header ending in a tab is rejected as having a blank column name.
- Plasmid tables get no duplicate / replicon-name / tier-vocabulary QA.

## Software boundaries

- Starts from manifests and tool outputs; it does **not** run read QC, assembly, annotation, AMR detection, typing or phylogeny.
- Not designed for patient-identifiable data on a public Space. No access control, encryption at rest or audit service.
- Single-process pandas; benchmarked to 20,000 isolates (60,000 AST rows) in ~14 s on the authoring machine, not tested beyond that.
- Label harmonisation folds case/whitespace only; taxonomy and drug-name synonym dictionaries belong upstream.
