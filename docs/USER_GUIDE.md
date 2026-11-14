# User guide

For installation and the three-line quick start see the [README](../README.md). This guide explains how to use the app and command line, and how to read the output.

## 1. Install and check

```bash
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r requirements.txt && pip install -e .
amr-surveil doctor
```

`amr-surveil doctor` checks your Python packages and runs a 24-isolate end-to-end test. Fix every `FAIL` before anything else.

**Common first-run problems**

| Symptom | Cause and fix |
|---|---|
| `ModuleNotFoundError: pyarrow.lib` in the web app | A broken binary wheel, usually after mixing conda and pip, or a Python newer than the wheels (e.g. 3.14). `pip uninstall -y pyarrow && pip install -U pyarrow`, or use Python 3.11–3.13: `conda create -n amr python=3.12`. |
| `ModuleNotFoundError: No module named 'rpds.rpds'` | A broken `rpds-py` binary wheel (Python 3.14 on Windows has hit this). The app itself no longer imports the altair chart stack, but other Streamlit elements can. `pip uninstall -y rpds-py && pip install -U rpds-py`, or use Python 3.11–3.13. `amr-surveil doctor` checks it. |
| `streamlit un app.py` / `run app.p` | Typos in the command: it is `streamlit run app.py`. |
| App reachable from other computers | The shipped config binds to `127.0.0.1` (this computer only). The Docker image binds `0.0.0.0` on purpose, as containers need it. |
| Deprecation warnings about `use_container_width` | Fixed in 1.1.0; upgrade. |

## 2. Start the app

```bash
streamlit run app.py
```

It opens on a synthetic demo (nothing real). The sidebar has the upload slots and analysis settings; the tabs show results.

## 3. Using your own data

1. Prepare the tables (see the [data dictionary](DATA_DICTIONARY.md)); `amr-surveil templates --out my_inputs` writes header-only templates.
2. Upload `samples.csv` (required) and any others. `ast.csv` gives resistance indicators; `amr_calls.csv` plus `screening.csv` and a cited `gene_drug_map.csv` enable genotype–phenotype analysis.
3. Choose the deduplication rule your protocol requires.
4. Tick the study-declaration boxes **only if true**; they feed the go/no-go gate.
5. Press **Run analysis**. In strict mode (default) ERROR-level data problems stop the run and are listed.

Equivalent command line: `amr-surveil run --samples samples.csv --ast ast.csv --config config.yaml --out results/`.

## 4. Using public data (Public data hub tab)

| Source | How | What you get |
|---|---|---|
| BV-BRC (live) | Organism, antibiotic(s), optional years/countries, record cap | Isolate-level AST joined to genome metadata |
| BV-BRC (file) | Upload an "AMR phenotypes" table exported from the website | Same, from your file |
| NCBI Pathogen Detection | Upload an Isolates Browser / metadata table with `AST_phenotypes` | Isolate-level AST (method/standard `unspecified`) |
| ECDC Surveillance Atlas | Upload an EARS-Net export containing counts | Aggregate counts (map and context stream) |

Each source shows its **terms, attribution and caveats** first. After retrieval, read the notes (what was dropped and why), then press **Analyse this dataset**; every tab now shows that dataset. **Download these inputs** saves canonical CSVs plus `public_data_provenance.json` (query URL, retrieval time, response hash).

Command-line equivalents:

```bash
amr-surveil public-data sources
amr-surveil public-data bvbrc --organism "Escherichia coli" --antibiotic ciprofloxacin --max-records 2000 --out public_ecoli
amr-surveil public-data import --source ncbi_pathogen --input isolates.tsv --out public_ncbi
amr-surveil run --samples public_ecoli/samples.csv --ast public_ecoli/ast.csv --out results_public
```

**Reading public data correctly.** Repository genomes are a convenience sample: outbreak and resistant isolates are over-represented, there is no patient ID (so no deduplication), and test methods vary between submitters. Use the output to explore and generate questions, not to estimate national prevalence or rank countries.

## 5. Tabs

* **AST indicators**: R/(S+I+R) with 95% Wilson intervals by group; `I` counts as not-resistant; small groups are flagged.
* **Antibiogram**: % susceptible = S/(S+I+R) (so `I` counts against susceptibility, unlike the resistance proportion), minimum 30 isolates by default, one row per method/breakpoint context.
* **Map**: choose organism, drug, sector and **one** test-method/breakpoint context. Counts are summed, never percentages averaged. Countries below the minimum (default 20) are listed with counts but not coloured. Unrecognised country names are listed, not guessed.
* **Comparability**: flags (method, breakpoint version, laboratory, volume, gaps). A flag means do not read a trend into the numbers yet.
* **Genotype–phenotype**: runs only if the 7-criterion gate says GO; otherwise it states which criteria failed.
* **AST re-check**: re-derives S/I/R from MIC/zone with *your* breakpoint table.
* **Model outputs**: predictions (e.g. EpiResNet) stay separate; agreement needs an attested independent held-out split.
* **Context streams**: aggregate counts, incidence, consumption; never joined to isolates.
* **Data quality / Inventory**: findings, exclusions, harmonised labels, missingness.
* **Run provenance**: manifest with hashes and versions; `amr-surveil verify results/` re-checks them.

## 6. Reproducibility and sharing

Keep the whole results folder. `amr-surveil verify RESULTS --check-inputs` proves nothing changed. Report the software version, deduplication rule, breakpoint standard/version, and for public data the retrieval date and source.

## 7. Governance

Do not upload patient-identifiable or restricted data to a public-hosted Space. The software does not authorise clinical or operational use.
