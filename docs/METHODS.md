# Methods

Every rule below is implemented in `amr_surveil/` and pinned by a test. Parameters that change a result are declared in
[`config/default.yaml`](../config/default.yaml), validated, and echoed into `run_manifest.json`.

## 1. Record eligibility (what enters an indicator)

An AST record contributes to an indicator only if all hold; otherwise it is written to `exclusion_log.csv` with a reason.

1. `result` is `S`, `I` or `R` (blank/invalid → `missing_or_invalid_result`).
2. Its `sample_id` is present **once** in the manifest (`ambiguous_sample_id`, `unmatched_sample_id`).
3. The AST organism label, if given, agrees with the manifest (`organism_mismatch`).
4. It survives the chosen deduplication rule (§3).
5. It is not an exact duplicate (collapsed once) or a *conflicting* duplicate — same isolate/drug/method/breakpoint with different
   categories (`conflicting_duplicate_ast`; neither is counted, because choosing one would be invention).

Case/whitespace variants of organism, drug, country, sector and method labels are folded to the most frequent spelling
(ties: alphabetical) **pooled across all tables**, so the same species cannot split in two because the AST file and the manifest
spell it differently; changes are listed in `label_harmonisation.csv`. Evidence-tier labels are compared case-insensitively.
Records without a `sample_id` or an antimicrobial are logged (`missing_sample_id`, `missing_antimicrobial`) and removed; rows that
the aggregate/incidence/consumption calculations cannot use are logged as `invalid_row_excluded_from_context`. Synonyms (*E. coli* vs *Escherichia coli*) are **not** merged.

## 2. Resistance proportion

`R / (S + I + R)` per group (year × country × sector × organism × drug × **AST method** × **breakpoint standard/version**).
`I` stays in the denominator and out of the numerator (EUCAST: *susceptible, increased exposure*). A group counts **isolates**
(one per `sample_id`), not records. Intervals are Wilson score (default 95%), validated against Newcombe (1998) reference values
in `tests/test_epidemiology.py`. Percentages are never averaged across groups. `small_denominator_warning` is `n < threshold`
(default 10) — a display flag, not a protocol rule.

## 3. Deduplication (off by default)

Rules: `none`; `first_per_patient_species_year`; `first_per_patient_species_specimen_year`. Surveillance protocols differ, and
the tool does not choose one for you. The rule is applied among isolates that contribute at least one eligible AST record, so an
untested first isolate cannot remove a tested repeat. The earliest collection date wins; on a tie the lowest `sample_id` is kept
and `date_tie` is set, because "first" is then arbitrary. A year- or month-only date cannot be ordered against a full date in the
same period, so any group containing one is also flagged `date_tie`. The rule applies to resistance indicators **and** to the
genotype–phenotype analysis (repeat isolates would otherwise inflate the 2×2 counts). `patient_id` is compared case-insensitively. Samples lacking patient, organism, valid date (or specimen) are **retained** and logged
`not_assignable`. A rule that is requested while no `patient_id` is populated raises an error rather than silently doing nothing.

## 4. Comparability flags

Per organism × drug × country × sector series (nothing is adjusted or smoothed):

| Flag | Trigger |
|---|---|
| `unknown_year_records` | records without a valid year cannot be placed in a trend |
| `single_year_only` (INFO) | fewer than two years |
| `small_denominator_years` | years with fewer than the small-n threshold of isolates |
| `missing_years` | gaps inside the observed range |
| `testing_volume_jump` | consecutive-year isolate counts differ by ≥ `testing_volume_jump_ratio` (default 2×) |
| `breakpoint_version_changes_across_years` | different standard/version sets in different years |
| `ast_method_changes_across_years` | different AST methods in different years |
| `testing_laboratory_set_changes` | different laboratory sets in different years |
| `specimen_mix_changes` | different specimen sets in different years |

## 5. AST category re-check

If a breakpoints table is supplied, each record is matched (case-insensitively) on organism, drug, standard and version, and the
category is re-derived from MIC (preferred) or zone diameter. Censoring is respected: `<=x`/`<x` can only establish S, `>x`/`>=x` only R,
and only when the bound itself is on the correct side of the breakpoint; otherwise the result is `indeterminate`, never a guess.
Breakpoints are in mg/L: a MIC is used only if `mic_unit` is blank (assumed mg/L) or numerically identical (`mg/L`, `µg/mL`);
any other unit gives `unit_not_mg_per_l`. Records with a missing standard or version never match a breakpoint row, and a blank AST
organism is taken from the manifest. A censored MIC that cannot be resolved stays `indeterminate` (no silent fall-back to a zone
measured by a different method).
Statuses: `concordant`, `discordant`, `indeterminate`, `no_breakpoint`, `no_measurement`, `unit_not_mg_per_l`. Breakpoints are **never built in** —
they change between versions and the table you supply is recorded by hash in the manifest.

## 6. Genotype–phenotype analysis

### 6.1 Go / no-go gate

All seven criteria must be met or concordance is not computed (the report lists the failures):

1. WGS screening and AST share isolates by exact `sample_id` (≥1 isolate with a single screening record and a valid S/I/R).
2. AST method recorded for every paired AST record.
3. Breakpoint standard **and** version recorded for every paired AST record.
4. Screening method, database and database version recorded for every paired isolate.
5. No organism conflicts between manifest and AST among paired isolates.
6. A gene→drug map was supplied and every row cites a source.
7. The analyst attests, in `study:`, that the sample-selection process is documented (software cannot verify this).

### 6.2 Concordance

*Genotype-positive* = the isolate (declared screened, one screening record) carries ≥1 gene mapped to that organism/drug by the cited map
(`exact`, or `prefix` where the map says so) among calls of an accepted evidence tier. *Genotype-negative* = screened and no such gene.
Unscreened isolates are excluded — a missing call row is never read as a negative.
**Uncertain genotypes:** a call that matches a mapped gene but was set aside (weak evidence tier, or an AMRFinderPlus `PARTIAL*` /
`INTERNAL_STOP*` hit) makes an isolate *neither* positive nor negative for that drug unless it also has an accepted mapped hit. Such
isolates are left out of the 2×2 table, logged (`uncertain_genotype_low_evidence_hit`) and counted in `n_excluded_uncertain_genotype`. Phenotype-positive = `R`; `I` and `S` are not-R
(`n_I_in_notR` is reported). Counts: TP gene+/R, FP gene+/not-R, FN gene−/R, TN gene−/not-R; sensitivity, specificity, PPV, NPV with
Wilson intervals. Strata with fewer than `concordance_min_isolates` report counts only. `genotype_phenotype_discordant.csv` lists
every FP/FN isolate so a person can investigate (unmapped mechanism, regulation/expression, breakpoint or method issue, mislabelled
sample, non-functional gene).

**Interpretation boundary:** descriptive concordance of a mapped-determinant screen with the recorded phenotype in the supplied
isolates. It is not diagnostic accuracy, a causal attribution, a prediction for new isolates, or a clinical statement, and it
inherits every selection bias of the sample.

## 7. Model predictions

Predictions are never written into the AST table. For each model × drug × evaluation split the tool joins unique predictions to unique
observed AST on exact ID + drug and reports the four counts. Several predictions for the same isolate/drug (EpiResNet predicts
*protein-antibiotic pairs*) are not collapsed into an isolate label; they are excluded and logged. Agreement percentages with Wilson
intervals are produced only if `study.predictions_independence_attested` is true **and** `evaluation_split` is `test`/`held_out`
**and** there are ≥10 pairs. AMR-Surveil cannot check training-set overlap; the attestation is the analyst's responsibility.

## 8. Separate streams

- **Aggregate resistance:** `resistant_n / tested_n` from source counts; duplicate keys are an error; never joined to isolates.
- **Incidence:** `resistant_infection_n / population_n × multiplier` (default per 100,000). Not derivable from isolate AST.
- **Consumption:** `ddd_total / population_n / period_days × 1000` = DDD per 1,000 inhabitants per day. A population indicator, not patient exposure.

## 9. Provenance

`run_manifest.json` stores software, Python and pandas versions, the effective (validated) configuration, and the SHA-256 of every
input and every output. Given the same inputs and configuration, outputs are byte-identical apart from the two files that contain a
timestamp (`analysis_summary.json`, `report.md`) — asserted by a test.
