# Integration guide

All adapters write the canonical CSV contracts in [`DATA_DICTIONARY.md`](DATA_DICTIONARY.md). Nothing is invented: when a source does not provide a field, the field stays `NA` and the audit says so.

## 0. Do the adapters work? (honest status)

| Adapter | Verified | Not verified |
|---|---|---|
| AMRFinderPlus TSV | v3-style and v4-style headers, header-only and empty files, duplicate samples, partial/internal-stop hits | A report from a v4 binary (aliases come from documentation) |
| FishAMR-Link | Feature-table, JSON-mapped export and SQLite warehouse paths, written from reading the upstream v17.0 source (commit `fa84114`) | A capture from a live FishAMR-Link session |
| EpiResNet | Held-out `test_predictions.csv`, metrics and manifest handling, written from the upstream v5 source (commit `4cc89fa`) | A capture from a real training run |

Fixtures are hand-written to match those layouts. **Please test your own real exports once**: run the `normalize-*` command, then `amr-surveil qa` on the result. If a column is not recognised the error names it; send the (header-only) file layout and the adapter can be corrected.

## AMR-UNBIND (deliberately not ingested)

AMR-UNBIND (v1.3.2) runs steered molecular dynamics for **one protein–ligand pair per simulation** and reports steering force, protocol work and force-slope descriptors (`ml_features.csv`, `smd_analysis.csv`, `run_metadata.json`). Its own documentation states that these are not binding free energies. AMR-Surveil's tables are **isolate-level** (sample, organism, drug, S/I/R, gene call). A simulation descriptor has no isolate, no organism-level denominator and no breakpoint, so joining it to AST records would imply a mechanistic or causal link that neither tool can support, and the AMR-Surveil rule that unlike evidence streams are never joined would be broken. There is no adapter. If you want to discuss simulations next to surveillance results, keep them as a separately cited table in your report, keyed by protein target and drug (for example, "TEM-1 + ampicillin"), and do not merge them into `samples.csv` or `ast_results.csv`.

## Using another tool of yours

Any tool can feed AMR-Surveil without a dedicated adapter: write a CSV with the columns of the target table in [`DATA_DICTIONARY.md`](DATA_DICTIONARY.md) (`amr-surveil templates --out DIR` creates header-only templates), or map its columns with a JSON mapping (`normalize-fishamr --kind amr --mapping` / `normalize-epiresnet --mapping`, examples in `config/`). A dedicated adapter is warranted when the tool's output has its own quirks (like AMRFinderPlus's header-only "screened, nothing found" files). To add one, share a small real output file and the tool's version.

## 1. AMRFinderPlus (recommended detection engine)

```bash
amr-surveil normalize-amrfinder \
  --input sampleA.tsv sampleB.tsv sampleC.tsv \
  --organism "Escherichia coli" --tool-version 3.12.8 --database-version 2024-01-31.1 \
  --output normalized/amr_calls.csv --screening-output normalized/screening.csv
```

- Header variants of v3.x and v4.x are mapped; element types default to `AMR` (override with `--element-types`).
- A report with a header and **no rows** means *screened, nothing found* → the sample is written to the screening table with zero calls. A **0-byte** file is an error, not a negative.
- The tool and database versions are not stored in AMRFinderPlus reports, so you must supply them; they are required by the gate.
- Pass `--sample-id` when the report has no `Name` column; two reports for the same sample are an error.

## 2. FishAMR-Link

What the adapters rely on (read from FishAMR-Link v17.0, commit `fa84114`): the app has no CSV export and its rows carry no sample ID; its fallback scanner reports fixed `confidence` constants (0.89 / 0.92, versus 1.0 for AMRFinderPlus) that are not statistical; the remote-BLAST fallback writes generic "Homology Match" hits; the warehouse seeds simulated reference isolates; `phenotypic_resistance` is free text.

```bash
# in-app feature table (needs an explicit sample ID because rows carry none)
amr-surveil normalize-fishamr --kind feature-table --input features.csv --sample-id ISO-001 \
  --organism "Escherichia coli" --output normalized/amr_calls.csv --plasmids-output normalized/plasmid_calls.csv
# SQLite warehouse (opened read-only; simulated baseline skipped; coordinates not exported)
amr-surveil normalize-fishamr --kind warehouse --input fishamr.db --output normalized/amr_calls.csv --samples-output normalized/samples.csv
```

Each row receives an **evidence tier** (`curated_database` for AMRFinderPlus-derived rows, `homology_search`, `heuristic_motif`, `unspecified`). The two lowest tiers are excluded from genotype analysis by default; fixed "confidence" values are ignored. Free-text phenotypes are **never** converted to AST. Run `amr-surveil normalize-fishamr --help` for the exact flag names of your version.

## 3. EpiResNet (v5; commit `4cc89fa`)

```bash
amr-surveil normalize-epiresnet --test-predictions --input test_predictions.csv \
  --final-metrics final_metrics.json --run-manifest run_manifest.json --checkpoint best.pt \
  --output normalized/model_predictions.csv
```

- `test_predictions.csv` is the **held-out test split** only, and rows are **protein–antibiotic pairs**; the adapter says so in every row (`prediction_unit`).
- Threshold, temperature, split strategy and training-data hash are copied from the metrics/manifest; the checkpoint is identified by SHA-256 without loading PyTorch.
- EpiResNet's training script does not export predicted MIC; `predicted_mic` stays empty unless your file has `predicted_log2_mic`.
- `isolate_id` becomes `sample_id`; make sure it lives in the same ID namespace as your WGS `sample_id`, or the predictions will be reported as unmatched.
- Agreement statistics need `study.predictions_independence_attested: true` **and** a held-out split. The tool cannot verify independence.

## 4. Your own lab files

Use `normalize-fishamr --kind amr|plasmid` or `normalize-epiresnet --mapping` (without `--test-predictions`) with a JSON column mapping (`config/*.example.json`) for any CSV export. Then:

```bash
amr-surveil run --samples samples.csv --amr-calls amr_calls.csv --screening screening.csv --ast ast.csv \
  --breakpoints breakpoints.csv --gene-drug-map gene_drug_map.csv --config config.yaml --out results/
```

Before the first real run read [`LIMITATIONS.md`](LIMITATIONS.md) and fill the `study:` block of the configuration honestly.
