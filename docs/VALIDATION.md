# Validation record

This page separates **what has been demonstrated** from **what would still be needed** before AMR-Surveil could be called
validated for any operational purpose. It is a software-verification record, not a biological validation.

## 1. Automated tests (verified)

`python -W error -m pytest -q tests` → **388 passed** (Python 3.13.18, pandas 3.0.5, warnings promoted to errors).

| Area | Representative evidence |
|---|---|
| Statistics | Wilson intervals vs published Newcombe (1998) reference values (81/263 → 25.53–36.62; 0/20 → 0–16.11; 20/20 → 83.89–100); small-n boundary |
| Concordance | Hand-computed 2×2 tables, `I` counted as not-R, screened-only negatives, suppression below the minimum stratum size |
| AST logic | Censored-MIC decisions incl. indeterminate cases; breakpoint matching; zone diameters |
| Validation rules | `tests/test_validation_matrix.py`: 30 QA codes, each with a triggering and a non-triggering input |
| Adapters | AMRFinderPlus v3/v4 headers, header-only (screened, zero calls), empty file, duplicate sample; FishAMR-Link feature table, warehouse (seeded baseline excluded, coordinates never exported); EpiResNet held-out predictions |
| Pipeline | End-to-end run, strict vs exploratory, deterministic outputs, manifest hashes of inputs *and* outputs |
| Contracts | JSON schemas regenerate identically; every synthetic table conforms; documentation, version and CLI help stay in sync (`test_docs_sync.py`) |
| UI | `app.py` executed end-to-end against stubbed `streamlit`/`plotly` (all tabs, upload, public-hub and failure paths); static AST test forbids bare conditional-expression statements (Streamlit magic) and deprecated arguments |
| Public data | Mocked-transport tests: HTTPS/host/credential refusal, size and gzip-bomb limits, cache tamper detection, BV-BRC/NCBI/ECDC parsers, country normalisation, map and antibiogram tables, public data through the full pipeline |

## 2. Mutation spot check (verified, limited)

Six faults were injected by hand into core logic and each was caught by at least one test after two weak spots were closed:

| Injected fault | Outcome |
|---|---|
| Wilson z constant altered | initially survived → reference-value test added → killed |
| Small-n boundary `<` vs `<=` | initially survived → boundary test added → killed |
| Concordance FP/FN swapped | killed |
| Censored-MIC operator direction flipped | killed |
| Evidence-tier filter disabled | killed |
| Go/no-go criterion bypassed | killed |

This is a spot check, not a full mutation-testing run.

## 3. Synthetic ground-truth benchmark (verified, synthetic only)

`benchmarks/run_benchmarks.py` generates data with a **planted** association (P(R | gene) = 0.92, P(R | no gene) = 0.06) and checks that the
pipeline recovers it, then records throughput. Reference run (`benchmarks/results/reference_20261111.json`, seed 20261107, 1 machine):

| Isolates | AST rows | Time (s) | Isolates/s | PPV (planted 0.92) | R without gene (planted 0.06) |
|---:|---:|---:|---:|---:|---:|
| 100 | 300 | 0.60 | 167 | 0.898 | 0.052 |
| 1,000 | 3,000 | 1.18 | 847 | 0.917 | 0.055 |
| 5,000 | 15,000 | 3.78 | 1,321 | 0.921 | 0.060 |
| 20,000 | 60,000 | 14.68 | 1,363 | 0.921 | 0.062 |

Convergence on the planted values shows the counting logic is right. It says nothing about real-world accuracy: the planted structure is
the only structure in the data.

## 4. Defects found in v0.1.0 (verified fixed)

v0.1.0 was probed adversarially; five classes of defect were reproduced **on the original code** and then fixed and pinned
(`tests/test_v010_regressions.py`): (A) blank AST result stopped strict runs as an ERROR; (B) censored MICs flagged invalid; (C) case variants of labels silently split strata;
(D) duplicate aggregate rows silently double-counted; (E) future / ambiguous dates accepted. Further fixes are listed in the [changelog](../CHANGELOG.md).

## 4b. Fuzzing, independent recomputation and independent review (verified)

- **Fuzzing** (`tests/test_robustness.py`, plus a larger 2,300-run campaign): inputs are damaged at random (blank/hostile cells, dropped columns,
  duplicated rows, empty tables, random deduplication rules, strict/exploratory). Every outcome must be a finished run or an
  `AMRSurveilError`. The campaign found one real crash (a plasmid row without `sample_id` under pandas 3), now fixed.
- **Independent implementation** (`tests/test_independent_oracle.py`): a second, stdlib-only (no pandas) computation of the resistance
  proportions, Wilson intervals and genotype–phenotype 2×2 counts from the raw CSVs agrees with the pipeline on five synthetic data sets.
- **Two adversarial review rounds** by an independent reader who had to demonstrate each defect with a running input: 17 defects were
  reproduced (deduplication, organism-label splitting, MIC units, breakpoint-key matching, low-evidence hits counted as negatives, silent
  drops, case-sensitive tiers, out-of-range dates, infinite values, false duplicate stops, mislabelled CI, adapter crashes). All are fixed with
  a regression test each ([changelog](../CHANGELOG.md)). Remaining, deliberately unfixed items are listed in [`LIMITATIONS.md`](LIMITATIONS.md).
  Reviews find what they find: this lowers, but cannot eliminate, the chance of further defects.
- **Clean-install simulation:** the Docker `COPY` set was copied into a clean directory, built into a non-editable install in a fresh
  environment and run from outside the source tree (`demo` and `run`). This is not a Docker build.

## 4c. Independent review of the 1.1.0 public-data code (verified)

A reviewer that had not seen the implementation tried to break the new modules and reproduced 11 defects (listed in the changelog), all fixed with regression tests. It also found no defect in the redirect/host/credential checks, gzip size cap, map pooling guards or `doctor`. The connectors were still not run against the live services.

## 5. Coverage

About 96% of lines (re-measured after the hardening round), measured with a stdlib `sys.settrace` script because `coverage.py` could not be installed offline. It is an approximation
(it was measured before the last validation-rule tests were added). CI runs `pytest-cov` with a floor; trust that number once CI is green.

## 6. What has **not** been validated

- Docker image build and container start (no daemon was available); CI contains build and health-check steps.
- Streamlit in a real browser for 1.1.0's new tabs; Plotly rendering of the choropleth.
- Public-data connectors against the live BV-BRC, NCBI and ECDC services.
- Python 3.10–3.12 and pandas 2.2.x (CI matrix only).
- Adapters on captures from live FishAMR-Link / AMRFinderPlus / EpiResNet runs; AMRFinderPlus v4 aliases are taken from documentation.
- Any real surveillance data set, any external benchmark, any laboratory-accredited workflow.
- Scientific review of the gene→drug logic by a microbiologist (the bundled map is a **synthetic fixture**, not a curated resource).

## 7. Plan to close the gaps

1. Green CI on the matrix (Python 3.10–3.13, pandas 2.2.3 and 3.x) including Docker build and health check.
2. Capture real outputs (AMRFinderPlus 3.x/4.x, FishAMR-Link feature table, EpiResNet test split) on public isolates and add them as fixtures.
3. Use a public paired WGS + AST collection (with its licence recorded) to run the gate and concordance end-to-end; document the gate outcome even if NO-GO.
4. Independent review of the methods by a microbiologist/epidemiologist.
5. Run `mutmut` or equivalent on the statistics, AST and concordance modules.
