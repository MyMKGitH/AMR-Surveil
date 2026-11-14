"""Streamlit interface for AMR-Surveil (demo and private/local exploratory runs).

All logic lives in ``amr_surveil.ui_support`` and the pipeline; this file only lays out widgets.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from streamlit.errors import StreamlitAPIException

from amr_surveil import __version__
from amr_surveil.antibiogram import cumulative_antibiogram
from amr_surveil.datasources import (
    SOURCES,
    FetchError,
    default_cache_dir,
    fetch_bvbrc,
    parse_bvbrc_export,
    parse_ecdc_atlas,
    parse_pathogen_table,
)
from amr_surveil.exceptions import AMRSurveilError
from amr_surveil.help_text import QUICK_START_MD
from amr_surveil.mapping import country_resistance
from amr_surveil.ui_support import (
    DEDUP_CHOICES,
    TIER_CHOICES,
    UPLOAD_SLOTS,
    analyse_public_dataset,
    build_config,
    context_labels,
    error_bars,
    filter_context,
    public_inputs_zip,
    read_result,
    read_table_upload,
    run_demo,
    run_from_files,
    zip_outputs,
)

st.set_page_config(page_title="AMR-Surveil", page_icon="🧬", layout="wide")


def _wide(widget, *args, **kwargs):
    """Call a Streamlit widget at full width on both old and new releases.

    Newer Streamlit deprecates ``use_container_width`` in favour of ``width="stretch"``; older
    releases reject the string or do not know the keyword. Trying the modern spelling first avoids
    the deprecation warning without pinning a minimum version.
    """
    try:
        return widget(*args, width="stretch", **kwargs)
    except (TypeError, StreamlitAPIException):
        return widget(*args, use_container_width=True, **kwargs)


def show_table(
    frame: pd.DataFrame, empty_message: str = "No rows to show.", *, success: bool = False
) -> None:
    """Render a table, or a message when it is empty (explicit ``if``; never a bare conditional)."""
    if frame.empty:
        notice = st.success if success else st.info
        notice(empty_message)
        return
    _wide(st.dataframe, frame)


st.markdown(
    """
    <style>
      .block-container {padding-top: 1.8rem; padding-bottom: 2.5rem; max-width: 1440px;}
      .hero {padding: 1.4rem 1.6rem; border: 1px solid rgba(120,140,160,.35); border-radius: 18px;
             background: linear-gradient(120deg, rgba(30,90,130,.10), rgba(70,120,120,.06));}
      .muted {color: #788797; font-size: .9rem;}
      div[data-testid="stMetric"] {border: 1px solid rgba(120,140,160,.28); padding: 12px; border-radius: 12px;}
    </style>
    <div class="hero">
      <h1>🧬 AMR-Surveil</h1>
      <p><b>Reproducible genomic and epidemiological surveillance of antimicrobial resistance</b></p>
      <p class="muted">Provenance · screened-sample genotypes · comparability flags · go/no-go gate · transparent limitations</p>
    </div>
    """,
    unsafe_allow_html=True,
)
st.caption(
    f"AMR-Surveil v{__version__} · Research/portfolio software · Not a clinical decision-support system"
)

with st.sidebar:
    st.header("Run workspace")
    st.warning(
        "Do not upload patient-identifiable, confidential, or otherwise restricted data to a public-hosted Space. "
        "Use synthetic or properly licensed public data here; process governed clinical data only in an authorised environment."
    )
    uploaded = {
        key: st.file_uploader(
            f"{label}{' (required)' if required else ''}", type=["csv"], key=f"up_{key}"
        )
        for key, label, required in UPLOAD_SLOTS
    }
    st.markdown("**Analysis settings**")
    strict = st.checkbox("Strict quality gate (stop on ERROR findings)", value=True)
    dedup = st.selectbox(
        "Deduplication rule", DEDUP_CHOICES, index=0, help="Must follow your surveillance protocol."
    )
    excluded = st.multiselect(
        "Evidence tiers excluded from genotype analysis",
        TIER_CHOICES,
        default=["heuristic_motif", "homology_search"],
    )
    st.markdown("**Study declaration**")
    study_name = st.text_input("Study name")
    selection_ok = st.checkbox("Sample-selection process is documented")
    selection_text = st.text_area("Describe how isolates were selected", height=80)
    preds_ok = st.checkbox(
        "Model predictions are independent of the model's training data (attested)"
    )
    run_clicked = _wide(
        st.button,
        "Run analysis",
        type="primary",
        disabled=uploaded["samples"] is None,
    )
    demo_clicked = _wide(st.button, "Reload synthetic demo")

if run_clicked:
    for key in ("amr_surveil_result_dir", "amr_surveil_summary"):
        st.session_state.pop(key, None)
    workspace = Path(tempfile.mkdtemp(prefix="amr_surveil_upload_"))
    paths: dict[str, Path] = {}
    for key, item in uploaded.items():
        if item is not None:
            target = workspace / f"{key}.csv"
            target.write_bytes(item.getvalue())
            paths[key] = target
    result_dir = workspace / "results"
    try:
        config = build_config(
            strict=strict,
            deduplication_rule=dedup,
            excluded_evidence_tiers=excluded,
            study_name=study_name,
            sample_selection_documented=selection_ok,
            selection_description=selection_text,
            predictions_independence_attested=preds_ok,
        )
        st.session_state["amr_surveil_summary"] = run_from_files(paths, config, result_dir)
        st.session_state["amr_surveil_result_dir"] = str(result_dir)
        st.success("Analysis completed. Review the QA table before interpreting any result.")
    except (AMRSurveilError, ValueError, OSError) as exc:
        st.error(str(exc))
        issues = read_result(result_dir, "audit_issues.csv")
        if not issues.empty:
            show_table(issues)

if demo_clicked or ("amr_surveil_result_dir" not in st.session_state and not run_clicked):
    summary, result_dir = run_demo(Path(tempfile.mkdtemp(prefix="amr_surveil_demo_")))
    st.session_state["amr_surveil_result_dir"] = str(result_dir)
    st.session_state["amr_surveil_summary"] = summary
    st.info(
        "Showing a deterministic synthetic example. These values are not real surveillance observations."
    )

if "amr_surveil_result_dir" not in st.session_state:
    st.info("Upload a sample manifest or reload the synthetic demo to begin.")
    st.stop()

result_dir = Path(st.session_state["amr_surveil_result_dir"])
summary = st.session_state["amr_surveil_summary"]

cols = st.columns(6)
for col, (label, value) in zip(
    cols,
    [
        ("Samples", summary.get("n_samples", 0)),
        ("Screened", summary.get("n_screened_samples", 0)),
        ("AMR calls used", summary.get("n_amr_calls_used", 0)),
        ("AST records", summary.get("n_ast_records", 0)),
        ("QA errors", summary.get("n_qa_errors", 0)),
        ("Excluded records", summary.get("n_excluded_records", 0)),
    ],
    strict=True,
):
    col.metric(label, value)
if summary.get("n_qa_errors") and not summary.get("strict_quality_gate"):
    st.error("Exploratory run with unresolved ERROR findings — do not communicate these results.")

tabs = st.tabs(
    [
        "AST indicators",
        "Antibiogram",
        "Map",
        "Comparability",
        "Genotype–phenotype",
        "AST re-check",
        "Model outputs",
        "Context streams",
        "Data quality",
        "Inventory",
        "Run provenance",
        "Public data hub",
        "Help",
    ]
)
(
    tab_ast,
    tab_abg,
    tab_map,
    tab_comp,
    tab_gp,
    tab_recheck,
    tab_pred,
    tab_ctx,
    tab_qa,
    tab_inv,
    tab_prov,
    tab_hub,
    tab_help,
) = tabs

with tab_ast:
    frame = read_result(result_dir, "resistance_by_group.csv")
    st.subheader("Isolate-level resistance proportion with 95% Wilson intervals")
    st.caption(
        "R / (S + I + R) per group. Different methods and breakpoints are never pooled. Not a country ranking."
    )
    if frame.empty:
        st.info("No eligible isolate-level AST records were provided.")
    else:
        options = sorted(frame["antimicrobial"].dropna().astype(str).unique())
        selected = st.multiselect("Antimicrobial", options, default=options[:1])
        shown = frame[frame["antimicrobial"].isin(selected)].copy()
        up, down = error_bars(shown)
        shown["_up"], shown["_down"] = up, down
        shown["year_label"] = shown["year"].astype(str)
        fig = px.scatter(
            shown,
            x="year_label",
            y="resistant_pct",
            color="method",
            symbol="source_sector",
            error_y="_up",
            error_y_minus="_down",
            hover_data=[
                "country",
                "organism",
                "n_resistant",
                "n_tested",
                "breakpoint_standard",
                "breakpoint_version",
                "small_denominator_warning",
            ],
            title="Descriptive resistance proportions (small-n groups carry wide intervals)",
        )
        fig.update_xaxes(type="category", title="Year")
        fig.update_yaxes(title="Resistant isolates (%)", range=[0, 100])
        _wide(st.plotly_chart, fig)
        show_table(shown.drop(columns=["_up", "_down", "year_label"]))
with tab_abg:
    indicators = read_result(result_dir, "resistance_by_group.csv")
    st.subheader("Cumulative antibiogram (% susceptible)")
    st.caption(
        "CLSI M39-style: S / (S + I + R), so intermediate counts as not susceptible (unlike the resistance "
        "proportion elsewhere). Rows below the minimum isolate count are flagged. Methods and breakpoint "
        "versions are never merged. Deduplication follows the rule chosen for the run."
    )
    if indicators.empty:
        st.info("No eligible isolate-level AST records were provided.")
    else:
        years = pd.to_numeric(indicators["year"], errors="coerce").dropna().astype(int)
        c1, c2, c3, c4 = st.columns(4)
        if years.empty:
            st.info("No collection years are recorded, so no year filter is applied.")
            year_lo = year_hi = None
        else:
            year_lo = int(
                c1.number_input("From year", int(years.min()), int(years.max()), int(years.min()))
            )
            year_hi = int(
                c2.number_input("To year", int(years.min()), int(years.max()), int(years.max()))
            )
            n_undated = int(pd.to_numeric(indicators["year"], errors="coerce").isna().sum())
            if n_undated:
                st.warning(
                    f"{n_undated} indicator group(s) have no collection year and are excluded from this table."
                )
        sector_options = ["(all)"] + sorted(
            indicators["source_sector"].dropna().astype(str).unique()
        )
        sector_pick = c3.selectbox("Sector", sector_options, key="abg_sector")
        min_n_abg = c4.number_input("Minimum isolates", 1, 1000, 30)
        abg = cumulative_antibiogram(
            indicators,
            year_from=year_lo,
            year_to=year_hi,
            source_sector=None if sector_pick == "(all)" else sector_pick,
            min_n=int(min_n_abg),
        )
        show_table(abg, "No isolates match this selection.")
        if not abg.empty:
            low_rows = int((~abg["meets_min_n"]).sum())
            if low_rows:
                st.warning(
                    f"{low_rows} row(s) have fewer than {int(min_n_abg)} isolates; read them as counts, not percentages."
                )
            st.download_button(
                "Download antibiogram CSV", abg.to_csv(index=False), "antibiogram.csv", "text/csv"
            )
with tab_map:
    st.subheader("Where is resistance concentrated?")
    st.warning(
        "Read with care: a map of isolates is a map of **who sampled and sequenced**, not of national "
        "prevalence. Countries below the minimum isolate count are not coloured. One organism, one drug, "
        "one test method and one sector per map - nothing is pooled across them."
    )
    stream = st.radio(
        "Data stream",
        ["Isolate-level (this run)", "Aggregate counts (imported, e.g. ECDC)"],
        horizontal=True,
    )
    if stream.startswith("Isolate"):
        source_frame = read_result(result_dir, "resistance_by_group.csv")
        res_col, test_col = "n_resistant", "n_tested"
    else:
        agg = st.session_state.get("hub_aggregate")
        if agg is None or agg.empty:
            agg = read_result(result_dir, "aggregate_surveillance_context.csv")
        source_frame = (
            agg.rename(columns={"resistant_n": "n_resistant", "tested_n": "n_tested"})
            if not agg.empty
            else agg
        )
        res_col, test_col = "n_resistant", "n_tested"
    if source_frame.empty:
        st.info(
            "No data in this stream yet. Run an analysis or import data in the Public data hub tab."
        )
    else:
        m1, m2, m3 = st.columns(3)
        org_pick = m1.selectbox(
            "Organism", sorted(source_frame["organism"].dropna().astype(str).unique())
        )
        sub = source_frame[source_frame["organism"] == org_pick]
        drug_pick = m2.selectbox(
            "Antimicrobial", sorted(sub["antimicrobial"].dropna().astype(str).unique())
        )
        sub = sub[sub["antimicrobial"] == drug_pick]
        sectors_here = (
            sorted(sub["source_sector"].dropna().astype(str).unique())
            if "source_sector" in sub
            else []
        )
        sector_map = m3.selectbox("Sector", sectors_here) if sectors_here else None
        contexts = context_labels(sub)
        context_pick = (
            st.selectbox("Test method | breakpoint standard | version", contexts)
            if contexts
            else None
        )
        sub = filter_context(sub, context_pick)
        yrs = pd.to_numeric(sub["year"], errors="coerce").dropna().astype(int)
        if yrs.empty:
            st.info("No year information is available for this selection.")
        else:
            y1, y2, y3 = st.columns(3)
            lo = y1.number_input(
                "From year", int(yrs.min()), int(yrs.max()), int(yrs.min()), key="map_lo"
            )
            hi = y2.number_input(
                "To year", int(yrs.min()), int(yrs.max()), int(yrs.max()), key="map_hi"
            )
            min_n_map = y3.number_input("Minimum isolates per country", 1, 1000, 20, key="map_min")
            try:
                table, unmapped = country_resistance(
                    sub,
                    organism=org_pick,
                    antimicrobial=drug_pick,
                    year_from=int(lo),
                    year_to=int(hi),
                    source_sector=sector_map,
                    min_n=int(min_n_map),
                    resistant_col=res_col,
                    tested_col=test_col,
                )
            except ValueError as exc:
                st.error(str(exc))
            else:
                drawable = table[table["iso3"].notna() & ~table["suppressed"].astype(bool)]
                if drawable.empty:
                    st.info("No country reaches the minimum isolate count for this selection.")
                else:
                    fig = px.choropleth(
                        drawable,
                        locations="iso3",
                        color="resistant_pct",
                        hover_name="country",
                        hover_data=["n_resistant", "n_tested", "wilson_low_pct", "wilson_high_pct"],
                        color_continuous_scale="Reds",
                        range_color=(0, 100),
                        labels={"resistant_pct": "Resistant (%)"},
                        title=f"{org_pick} - {drug_pick}: resistant isolates (%), {int(lo)}-{int(hi)}",
                    )
                    fig.update_geos(showcountries=True, showframe=False)
                    _wide(st.plotly_chart, fig)
                show_table(table, "No countries to list.")
                if unmapped:
                    st.warning(f"Country names not recognised, so not drawn: {', '.join(unmapped)}")
                st.caption(
                    "Wilson 95% intervals are in the table. Suppressed countries show counts only."
                )
with tab_comp:
    frame = read_result(result_dir, "comparability_flags.csv")
    st.subheader("Threats to comparability across years")
    st.caption(
        "A flag means: do not read a trend into these numbers until the cause is understood. Nothing is adjusted."
    )
    show_table(frame, "No flags raised.", success=True)
with tab_gp:
    readiness = json.loads((result_dir / "readiness_assessment.json").read_text(encoding="utf-8"))
    st.subheader(f"Go / no-go decision: {readiness['decision']}")
    show_table(pd.DataFrame(readiness["criteria"]))
    conc = read_result(result_dir, "genotype_phenotype_concordance.csv")
    if readiness["decision"] != "GO":
        st.warning("Concordance was not calculated. Resolve the unmet criteria first.")
    elif not conc.empty:
        st.caption(
            "Descriptive concordance among screened isolates; not diagnostic accuracy. Small strata show counts only."
        )
        show_table(conc)
        st.subheader("Discordant isolate–drug records")
        show_table(read_result(result_dir, "genotype_phenotype_discordant.csv"))
    st.subheader("Exact-ID co-occurrence inventory")
    show_table(read_result(result_dir, "matched_isolates.csv"))
    st.subheader("Plasmid / replicon evidence")
    show_table(read_result(result_dir, "plasmid_calls_linked.csv"))
with tab_recheck:
    frame = read_result(result_dir, "ast_category_recheck.csv")
    st.subheader("Recorded categories re-derived from MIC/zone and supplied breakpoints")
    if frame.empty:
        st.info("Supply a breakpoints table to enable this check.")
    else:
        counts = frame["status"].value_counts().rename_axis("status").reset_index(name="rows")
        fig = px.bar(counts, x="status", y="rows", title="Recorded vs re-derived category")
        fig.update_xaxes(type="category")
        _wide(st.plotly_chart, fig)
        show_table(frame[frame["status"] != "concordant"])
with tab_pred:
    st.subheader("Model outputs vs observed AST")
    st.caption(
        "Predictions stay separate from AST. Accuracy-type values need an attested independent, held-out split."
    )
    show_table(read_result(result_dir, "prediction_vs_observed.csv"))
    show_table(read_result(result_dir, "model_prediction_context.csv"))
with tab_ctx:
    for title, name in (
        (
            "Aggregate resistance (source counts; not isolate-linked)",
            "aggregate_surveillance_context.csv",
        ),
        ("Population incidence (cases per stated denominator)", "incidence_context.csv"),
        (
            "Antimicrobial consumption (DDD per 1,000 inhabitants per day)",
            "consumption_context.csv",
        ),
    ):
        st.subheader(title)
        show_table(read_result(result_dir, name))
with tab_qa:
    frame = read_result(result_dir, "audit_issues.csv")
    st.subheader("Quality findings")
    show_table(
        frame,
        "No configured QA findings were raised. This does not certify the source data.",
        success=True,
    )
    st.subheader("Exclusion log")
    show_table(read_result(result_dir, "exclusion_log.csv"))
    st.subheader("Label variants merged")
    show_table(read_result(result_dir, "label_harmonisation.csv"))
with tab_inv:
    st.subheader("Data inventory")
    show_table(read_result(result_dir, "data_inventory.csv"))
    st.subheader("AST missingness")
    show_table(read_result(result_dir, "ast_missingness.csv"))
with tab_prov:
    manifest = json.loads((result_dir / "run_manifest.json").read_text(encoding="utf-8"))
    st.subheader("Run manifest")
    st.json(manifest)
    st.markdown((result_dir / "report.md").read_text(encoding="utf-8"))
    st.download_button(
        "Download results ZIP",
        data=zip_outputs(result_dir),
        file_name="amr_surveil_results.zip",
        mime="application/zip",
    )

with tab_hub:
    st.subheader("Public data hub")
    st.info(
        "Bring openly available data into the same pipeline. Data are fetched to this session on your "
        "request, never bundled with the software. Use public data only - not patient-identifiable data."
    )
    label_to_id = {
        "BV-BRC: query live (isolate-level AST + genome metadata)": "bvbrc",
        "BV-BRC: upload an AMR phenotype export": "bvbrc_file",
        "NCBI Pathogen Detection: upload an Isolates Browser / metadata table": "ncbi_pathogen",
        "ECDC Surveillance Atlas: upload an EARS-Net export (aggregate)": "ecdc_atlas",
    }
    choice = label_to_id[st.selectbox("Source", list(label_to_id))]
    info = SOURCES["bvbrc" if choice == "bvbrc_file" else choice]
    with st.expander("Terms, attribution and caveats for this source", expanded=True):
        st.markdown(f"**{info.name}** - {info.provider}  \n{info.scope}")
        st.markdown(
            f"**Terms** ({info.terms_checked}): {info.terms}  \n[Terms page]({info.terms_url})"
        )
        st.markdown(f"**Attribution:** {info.attribution}")
        for caveat in info.caveats:
            st.markdown(f"- {caveat}")
    cache_dir = default_cache_dir()
    dataset = None
    try:
        if choice == "bvbrc":
            h1, h2 = st.columns(2)
            organism_q = h1.text_input("Organism (species)", "Escherichia coli")
            drugs_q = h2.text_input("Antibiotic(s), comma-separated", "ciprofloxacin")
            h3, h4, h5 = st.columns(3)
            limit_q = h3.number_input("Max AST records", 100, 25000, 2000, step=100)
            yfrom = h4.number_input("Year from (0 = any)", 0, 2100, 0)
            yto = h5.number_input("Year to (0 = any)", 0, 2100, 0)
            countries_q = st.text_input("Only these countries (comma-separated, optional)")
            refresh_q = st.checkbox("Ignore cached copy and download again")
            if st.button("Retrieve from BV-BRC"):
                with st.spinner("Querying BV-BRC..."):
                    dataset = fetch_bvbrc(
                        organism_q,
                        [d for d in drugs_q.split(",") if d.strip()],
                        cache_dir=cache_dir,
                        max_records=int(limit_q),
                        year_from=int(yfrom) or None,
                        year_to=int(yto) or None,
                        countries=[c for c in countries_q.split(",") if c.strip()] or None,
                        refresh=refresh_q,
                    )
        else:
            upload = st.file_uploader(
                "Table (CSV or TSV)", type=["csv", "tsv", "txt"], key=f"hub_{choice}"
            )
            with_geno = False
            if choice == "ncbi_pathogen":
                with_geno = st.checkbox(
                    "Also import AMR genotypes (exploratory; empty cells are not treated as screened-negative)"
                )
            if upload is not None and st.button("Read file"):
                table_in = read_table_upload(upload.getvalue(), upload.name)
                if choice == "bvbrc_file":
                    dataset = parse_bvbrc_export(table_in)
                elif choice == "ncbi_pathogen":
                    dataset = parse_pathogen_table(table_in, include_genotypes=with_geno)
                else:
                    dataset = parse_ecdc_atlas(table_in)
    except (FetchError, AMRSurveilError, ValueError, OSError) as exc:
        st.error(str(exc))
    if dataset is not None:
        st.session_state["hub_dataset"] = dataset
        st.session_state["hub_aggregate"] = (
            dataset.aggregate if dataset.aggregate is not None else None
        )
    held = st.session_state.get("hub_dataset")
    if held is not None:
        st.success(
            f"Retrieved: {len(held.samples)} isolates, {len(held.ast)} AST records"
            + (
                f", {len(held.aggregate)} aggregate rows"
                if held.aggregate is not None and not held.aggregate.empty
                else ""
            )
        )
        for note in held.notes:
            st.warning(note)
        show_table(held.samples.head(50), "No isolate rows.")
        show_table(held.ast.head(50), "No AST rows.")
        if held.aggregate is not None:
            show_table(held.aggregate.head(50), "No aggregate rows.")
            st.caption(
                "Aggregate rows feed the Map tab (aggregate stream) and are never joined to isolates."
            )
        st.json(held.provenance, expanded=False)
        if not held.samples.empty and not held.ast.empty:
            if _wide(st.button, "Analyse this dataset", type="primary"):
                workdir = Path(tempfile.mkdtemp(prefix="amr_surveil_public_"))
                try:
                    cfg = build_config(
                        strict=False,
                        deduplication_rule="none",
                        excluded_evidence_tiers=excluded,
                        study_name=study_name
                        or f"Public data: {held.provenance.get('source_name', held.source_id)}",
                        sample_selection_documented=False,
                        selection_description="Convenience sample from a public repository; no sampling frame.",
                    )
                    summary_new, out_new = analyse_public_dataset(held, workdir, cfg)
                except (AMRSurveilError, ValueError, OSError) as exc:
                    st.error(str(exc))
                else:
                    st.session_state["amr_surveil_result_dir"] = str(out_new)
                    st.session_state["amr_surveil_summary"] = summary_new
                    st.rerun()
            st.download_button(
                "Download these inputs as canonical CSVs (with provenance)",
                public_inputs_zip(held, Path(tempfile.mkdtemp(prefix="amr_surveil_pubzip_"))),
                "public_dataset_inputs.zip",
                "application/zip",
            )
with tab_help:
    st.markdown(QUICK_START_MD)

st.divider()
st.caption(
    "AMR-Surveil is a research/portfolio workbench. It does not authorise diagnosis, treatment, prescribing, "
    "clinical reporting, or operational surveillance use without appropriate supervision, governance, and validation."
)
