# Roadmap

## Delivered
- **0.1.0** – descriptive integration layer: canonical contracts, QA audit, Wilson indicators, separate streams, Streamlit app, Docker.
- **0.2.0** – auditable genotype–phenotype workbench: screening semantics, evidence tiers, go/no-go gate, concordance with CIs, comparability flags, breakpoint re-check, real-format adapters (AMRFinderPlus, FishAMR-Link, EpiResNet v5), full provenance manifest, ground-truth benchmark.
- **1.0.0** – stable contract ([`STABILITY.md`](STABILITY.md)), `verify` and `templates` commands, fuzzing, an independent recomputation, two review rounds (17 defects fixed), 322 tests.
- **1.1.0** – first-run fixes from a real local run, `doctor`, public data hub (BV-BRC, NCBI Pathogen Detection, ECDC Atlas), map, cumulative antibiogram, user guide, related-work review.

## Next (ordered by evidence value)
0. **Run the public-data connectors against the live services** and fix any field-name drift; add captured responses as fixtures.
1. **Green CI + Docker** across Python 3.10–3.13 and pandas 2.2/3.x.
2. **Real fixtures**: captured AMRFinderPlus 3.x/4.x, FishAMR-Link and EpiResNet outputs on public isolates.
3. **Public paired WGS + AST case study** (licence recorded), reporting the gate decision and concordance, including failures.
4. **Independent review** of methods by a microbiologist / epidemiologist.
5. **Typing and phylogeny import** (MLST/cgMLST/SNP distance tables) with explicit interpretation limits.
6. **Raw-read workflow** (Nextflow/Snakemake) for one species and one toolchain, once chosen and tested.
7. **WHONET file import, MDR / co-resistance profiles** (needs a stated MDR definition), WHO GLASS import.
8. **R/Quarto companion report** and a SQL-backed store, if a use case needs them.

Milestones describe software; they are not claims of operational surveillance experience.
