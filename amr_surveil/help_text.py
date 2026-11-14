"""Built-in quick-start text shown in the app's Help tab (the full guide is docs/USER_GUIDE.md)."""

QUICK_START_MD = """
### 60-second tour
1. **Look first.** The app opens on a deterministic *synthetic* demo. Nothing in it is real data.
2. **Use your own data.** Upload `samples.csv` (required) plus any of the other tables in the sidebar,
   tick the study declaration boxes only if they are true, and press **Run analysis**.
3. **Or start from public data.** Open **Public data hub**, retrieve a dataset, read its terms and
   caveats, then press **Analyse this dataset**. Results appear in every other tab.
4. **Read the QA tab before any result.** A clean audit is not certification of the data.

### What each tab answers
| Tab | Question |
|---|---|
| AST indicators | What share of isolates is resistant, with what uncertainty, per drug/year/method? |
| Antibiogram | What % of isolates are susceptible (CLSI M39 style), and are there enough isolates? |
| Map | Where is resistance concentrated (one organism, one drug, one test method at a time)? |
| Comparability | Is a trend trustworthy, or did methods/labs/volumes change? |
| Genotype–phenotype | Do genes and AST agree? Only runs when the go/no-go gate says GO. |
| AST re-check | Do recorded S/I/R calls match MIC/zone + your breakpoint table? |
| Model outputs | How do predictions (e.g. EpiResNet) compare with observed AST? Kept separate. |
| Context streams | Aggregate counts, incidence and consumption; never joined to isolates. |
| Data quality / Inventory | What was excluded, merged or missing, and why? |
| Run provenance | Hashes, versions and settings needed to reproduce the run. |

### Rules the software enforces for you
* **I (intermediate) is counted as not-resistant**; resistance is R / (S + I + R).
* **Different methods and breakpoint versions are never pooled** - on a map, an antibiogram or a trend.
* **A missing gene is not "susceptible"**: negatives exist only for samples declared as screened.
* **Small groups are flagged or suppressed**, not shaded.
* **Public repository data are convenience samples**: outbreak and resistant isolates are
  over-represented. Proportions are not national prevalence and countries must not be ranked from them.

### If something fails
Run `amr-surveil doctor` in the same terminal/environment. It checks your Python packages and runs a
tiny end-to-end test. Missing `pyarrow` or an interpreter newer than the packages' wheels (for example
Python 3.14) are the usual first-run causes.

*Research software. Not a clinical decision-support system.*
"""
