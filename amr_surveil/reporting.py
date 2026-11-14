"""Human-readable Markdown report generated entirely from run outputs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def _cell(value: object) -> str:
    if value is None or (not isinstance(value, (list, dict)) and pd.isna(value)):
        return "—"
    return str(value).replace("|", "\\|").replace("\n", " ")


def _table(
    frame: pd.DataFrame, columns: list[str], headers: list[str], limit: int = 25
) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for row in frame.head(limit)[columns].itertuples(index=False):
        lines.append("| " + " | ".join(_cell(v) for v in row) + " |")
    if len(frame) > limit:
        lines += ["", f"*Showing {limit} of {len(frame)} rows; see the CSV output for all rows.*"]
    return lines


def _pct(value: object) -> str:
    return "—" if pd.isna(value) else f"{float(value):.1f}%"


def build_report(
    *,
    summary: dict[str, object],
    indicators: pd.DataFrame,
    flags: pd.DataFrame,
    issues: pd.DataFrame,
    readiness: dict[str, object],
    concordance: pd.DataFrame,
    discordant: pd.DataFrame,
    recheck: pd.DataFrame,
    missingness: pd.DataFrame,
    exclusions: pd.DataFrame,
    harmonisation: pd.DataFrame,
    pred_eval: pd.DataFrame,
    aggregates: pd.DataFrame,
    incidence: pd.DataFrame,
    consumption: pd.DataFrame,
    output_path: str | Path,
) -> None:
    """Write the report. Numbers come from the supplied frames and summary, never from prose."""
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    s = summary
    study = s.get("study") or {}
    L: list[str] = [
        "# AMR-Surveil analysis report",
        "",
        f"**Software version:** `{s.get('software_version', 'unknown')}`  ",
        f"**Run timestamp (UTC):** `{s.get('created_at_utc', 'unknown')}`  ",
        f"**Reference date for date checks:** `{s.get('reference_date', 'unknown')}`  ",
        f"**Deduplication rule:** `{s.get('deduplication_rule', 'none')}`  ",
        f"**Quality gate:** `{'strict' if s.get('strict_quality_gate') else 'exploratory (errors do not stop the run)'}`",
        "",
        "> **Scope:** research and portfolio analysis of supplied data. Not a clinical decision-support system, "
        "a diagnostic report, a national surveillance submission, or evidence of operational validation.",
        "",
    ]
    if not s.get("strict_quality_gate") and s.get("n_qa_errors"):
        L += [
            f"> ⚠️ **Exploratory run with {s['n_qa_errors']} unresolved ERROR finding(s).** "
            "Do not communicate these results until the errors are resolved.",
            "",
        ]
    if study.get("name") or study.get("question"):
        L += [
            "## Study definition",
            "",
            f"- **Name:** {_cell(study.get('name') or None)}",
            f"- **Question:** {_cell(study.get('question') or None)}",
            f"- **Inclusion criteria:** {_cell(study.get('inclusion_criteria') or None)}",
            f"- **Exclusion criteria:** {_cell(study.get('exclusion_criteria') or None)}",
            f"- **Sample-selection description:** {_cell(study.get('selection_description') or None)}",
            "",
        ]

    L += [
        "## Dataset inventory",
        "",
        "| Stream | Rows | Interpretation |",
        "|---|---:|---|",
        f"| Isolate/sample manifest | {s.get('n_samples', 0)} | Isolate-level metadata |",
        f"| Screened samples (genotyping declared) | {s.get('n_screened_samples', 0)} | Only these can have a *negative* genotype |",
        f"| AMR determinant calls (supplied / used) | {s.get('n_amr_calls', 0)} / {s.get('n_amr_calls_used', 0)} | Positive calls; evidence tiers excluded: `{', '.join(s.get('excluded_evidence_tiers', [])) or 'none'}` |",
        f"| Phenotypic AST records | {s.get('n_ast_records', 0)} | Observed categories supplied by the user |",
        f"| Aggregate surveillance rows | {s.get('n_aggregate_rows', 0)} | Source counts; not isolate-linked |",
        f"| Population incidence rows | {s.get('n_incidence_rows', 0)} | Cases per stated denominator; not isolate-linked |",
        f"| Antimicrobial consumption rows | {s.get('n_consumption_rows', 0)} | DDD-based indicator; not patient-level use |",
        f"| Model predictions | {s.get('n_predictions', 0)} | Research outputs; not observed AST |",
        f"| Plasmid/replicon rows | {s.get('n_plasmid_calls', 0)} | Upstream evidence; not proof of transfer |",
        f"| Isolates with AMR calls and AST (exact ID) | {s.get('n_matched_isolates', 0)} | Co-occurrence inventory only |",
        f"| Records excluded with a recorded reason | {s.get('n_excluded_records', 0)} | See `exclusion_log.csv` |",
        "",
    ]

    L += [
        "## Data-quality review",
        "",
        f"Errors **{s.get('n_qa_errors', 0)}** · Warnings **{s.get('n_qa_warnings', 0)}** · Info **{s.get('n_qa_info', 0)}**",
        "",
    ]
    if not issues.empty:
        L += _table(
            issues,
            ["severity", "table", "code", "row_count", "message"],
            ["Severity", "Table", "Code", "Rows", "Finding"],
            limit=40,
        )
        L.append("")
    if not harmonisation.empty:
        L += [
            "### Label variants merged for analysis",
            "",
            "Case/whitespace variants were folded onto the most frequent spelling. Synonyms (e.g. *E. coli* vs *Escherichia coli*) are **not** merged.",
            "",
        ]
        L += _table(
            harmonisation,
            ["table", "column", "canonical", "merged_variants"],
            ["Table", "Column", "Canonical", "Variants merged"],
            limit=15,
        )
        L.append("")
    if not exclusions.empty:
        L += ["### Exclusion summary", ""]
        L += _table(
            exclusions,
            ["stage", "table", "reason", "n_records"],
            ["Stage", "Table", "Reason", "Records"],
            limit=30,
        )
        L.append("")

    L += ["## Isolate-level AST resistance indicators", ""]
    if indicators.empty:
        L.append("No eligible AST records were available for indicator calculation.")
    else:
        L += [
            "Resistance proportion is **R / (S + I + R)** per group. `I` is kept in the denominator and not in the numerator "
            "(EUCAST: *susceptible, increased exposure*). Wilson score intervals. The small-n flag is a display warning, not a protocol rule.",
            "",
        ]
        disp = indicators.sort_values("n_tested", ascending=False).copy()
        disp["ci"] = [
            f"{r.resistant_pct:.1f}% ({r.wilson_low_pct:.1f}–{r.wilson_high_pct:.1f})"
            for r in disp.itertuples()
        ]
        disp["bp"] = (
            disp["breakpoint_standard"].astype(str) + " " + disp["breakpoint_version"].astype(str)
        )
        disp["rt"] = disp["n_resistant"].astype(str) + "/" + disp["n_tested"].astype(str)
        disp["caution"] = disp["small_denominator_warning"].map({True: "Small n", False: "—"})
        L += _table(
            disp,
            [
                "year",
                "country",
                "source_sector",
                "organism",
                "antimicrobial",
                "bp",
                "rt",
                "ci",
                "caution",
            ],
            [
                "Year",
                "Country",
                "Sector",
                "Organism",
                "Antimicrobial",
                "Breakpoint",
                "R / tested",
                f"Resistance ({round(100 * float(summary.get('confidence_level', 0.95)))}% CI)",
                "Caution",
            ],
        )
    L.append("")

    L += ["## Comparability flags", ""]
    if flags.empty:
        L.append(
            "No comparability threats were detected in the supplied data (this is not a guarantee of comparability)."
        )
    else:
        L.append(
            "Flags mean *do not read a trend into these numbers until the cause is understood*. Nothing is adjusted or smoothed."
        )
        L.append("")
        L += _table(
            flags,
            ["severity", "organism", "antimicrobial", "country", "source_sector", "flag", "detail"],
            ["Severity", "Organism", "Antimicrobial", "Country", "Sector", "Flag", "Detail"],
            limit=30,
        )
    L.append("")

    L += ["## AST category re-check against supplied breakpoints", ""]
    if recheck.empty:
        L.append(
            "No breakpoint table was supplied, or no AST record could be matched to a breakpoint. The recorded categories were not re-derived."
        )
    else:
        counts = recheck["status"].value_counts()
        L.append(
            "Categories re-derived from MIC/zone and the **supplied** breakpoint rows (never built in): "
            + ", ".join(f"{k} **{v}**" for k, v in counts.items())
            + "."
        )
        bad = recheck[recheck["status"] == "discordant"]
        if not bad.empty:
            L += ["", "Discordant records:", ""]
            L += _table(
                bad,
                [
                    "sample_id",
                    "antimicrobial",
                    "recorded_result",
                    "derived_result",
                    "mic",
                    "zone_mm",
                    "breakpoint_standard",
                    "breakpoint_version",
                ],
                [
                    "Sample",
                    "Antimicrobial",
                    "Recorded",
                    "Derived",
                    "MIC",
                    "Zone",
                    "Standard",
                    "Version",
                ],
                limit=20,
            )
    L.append("")

    L += [
        "## Genotype–phenotype analysis",
        "",
        f"**Go / no-go decision: `{readiness.get('decision')}`** — {readiness.get('n_paired_isolates', 0)} isolate(s) with a single screening record and a valid AST result.",
        "",
        "| Criterion | Met | Detail |",
        "|---|:---:|---|",
    ]
    for c in readiness.get("criteria", []):  # type: ignore[union-attr]
        L.append(
            f"| {_cell(c['description'])} | {'✅' if c['met'] else '❌'} | {_cell(c['detail'])} |"
        )
    L.append("")
    if readiness.get("decision") != "GO":
        L += [
            "Concordance was **not** calculated because the go/no-go conditions are not all met. Resolve the unmet criteria above first.",
            "",
        ]
    elif concordance.empty:
        L += [
            "The gate passed but no stratum could be formed (check the gene–drug map and organism/drug names).",
            "",
        ]
    else:
        L += [
            "Descriptive concordance of a mapped-determinant screen with recorded phenotype among **screened** isolates. "
            "Not diagnostic accuracy or causal attribution. Strata below the minimum isolate count show counts only.",
            "",
        ]
        view = concordance.copy()
        for col in ("sensitivity_pct", "specificity_pct", "ppv_pct", "npv_pct"):
            view[col] = view[col].map(_pct)
        L += _table(
            view,
            [
                "organism",
                "antimicrobial",
                "method",
                "breakpoint_standard",
                "breakpoint_version",
                "n_isolates",
                "n_gene_pos_R",
                "n_gene_pos_notR",
                "n_gene_neg_R",
                "n_gene_neg_notR",
                "sensitivity_pct",
                "specificity_pct",
                "ppv_pct",
                "npv_pct",
            ],
            [
                "Organism",
                "Drug",
                "Method",
                "Std",
                "Ver",
                "n",
                "G+/R",
                "G+/notR",
                "G−/R",
                "G−/notR",
                "Sens",
                "Spec",
                "PPV",
                "NPV",
            ],
        )
        if not discordant.empty:
            L += [
                "",
                f"{len(discordant)} discordant isolate–drug record(s) are listed in `genotype_phenotype_discordant.csv` for review "
                "(possible causes: unmapped mechanism, expression/regulation, breakpoint or method issue, mislabelled sample, gene not functional).",
            ]
    L.append("")

    L += ["## Model predictions versus observed AST", ""]
    if pred_eval.empty:
        L.append(
            "No predictions overlapped observed AST on exact sample ID and drug, or none were supplied."
        )
    else:
        L.append(
            "Predictions stay a separate stream. AMR-Surveil cannot verify that these isolates were absent from model training, "
            "so accuracy-type values appear only when the analyst attests independence **and** the rows are from a held-out split."
        )
        L.append("")
        view = pred_eval.copy()
        view["concordance"] = [_pct(v) for v in view["concordance_pct"]]
        L += _table(
            view,
            [
                "model_name",
                "antimicrobial",
                "evaluation_split",
                "n_pairs",
                "n_pred_R_obs_R",
                "n_pred_R_obs_notR",
                "n_pred_S_obs_R",
                "n_pred_S_obs_notR",
                "concordance",
                "metrics_status",
            ],
            [
                "Model",
                "Drug",
                "Split",
                "Pairs",
                "P-R/O-R",
                "P-R/O-notR",
                "P-S/O-R",
                "P-S/O-notR",
                "Agreement",
                "Status",
            ],
        )
    L.append("")

    L += ["## Missingness", ""]
    if missingness.empty:
        L.append("No AST records.")
    else:
        worst = missingness[missingness["n_missing_or_invalid_result"] > 0].sort_values(
            "pct_missing_or_invalid", ascending=False
        )
        if worst.empty:
            L.append("Every AST record has a valid S/I/R result.")
        else:
            L += _table(
                worst,
                [
                    "country",
                    "organism",
                    "antimicrobial",
                    "year",
                    "testing_lab",
                    "n_records",
                    "n_missing_or_invalid_result",
                    "pct_missing_or_invalid",
                ],
                [
                    "Country",
                    "Organism",
                    "Antimicrobial",
                    "Year",
                    "Lab",
                    "Records",
                    "Missing/invalid",
                    "%",
                ],
                limit=15,
            )
    L.append("")

    L += ["## Separate surveillance streams", ""]
    L.append(
        f"- **Aggregate resistance:** {len(aggregates)} row(s); percentages recomputed from source counts; never joined to isolates."
        if len(aggregates)
        else "- **Aggregate resistance:** none supplied."
    )
    L.append(
        f"- **Population incidence:** {len(incidence)} row(s); source-defined cases per stated population; resistance proportion ≠ incidence."
        if len(incidence)
        else "- **Population incidence:** none supplied."
    )
    L.append(
        f"- **Antimicrobial consumption:** {len(consumption)} row(s); DDD per 1,000 inhabitants per day; not patient-level exposure."
        if len(consumption)
        else "- **Antimicrobial consumption:** none supplied."
    )
    L += [
        "",
        "## Interpretation checklist",
        "",
        "- Confirm licence, accession provenance, retrieval date and the sample-selection process.",
        "- Verify metadata and AST belong to the same isolate; record method and breakpoint standard/version.",
        "- Review laboratory participation, specimen mix, missingness, deduplication, panel changes and breakpoint revisions before comparing years.",
        "- Do not compare incompatible human, animal, food, aquatic or environmental denominators as one sampling frame.",
        "- A negative gene call is not susceptibility; genomic similarity is not transmission.",
        "- Review all QA findings and `run_manifest.json` (it hashes every input and output) before communicating results.",
        "",
        "---",
        "Generated by AMR-Surveil. This report describes the supplied data and configuration; it does not certify data quality, clinical validity, laboratory performance or operational readiness.",
        "",
    ]
    target.write_text("\n".join(L), encoding="utf-8")
