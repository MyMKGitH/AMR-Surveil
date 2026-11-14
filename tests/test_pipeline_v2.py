import hashlib
import json

import pandas as pd
import pytest
import yaml

from amr_surveil import __version__
from amr_surveil.cli import DEMO_CONFIG
from amr_surveil.exceptions import ConfigurationError, QualityGateError
from amr_surveil.pipeline import OUTPUT_FILES, run_pipeline
from amr_surveil.synthetic import GENE_FOR_DRUG, PLANTED, generate_synthetic_dataset
from tests.conftest import pipeline_kwargs


def test_every_documented_output_exists_and_is_hashed(demo):
    out = demo["out"]
    manifest = json.loads((out / "run_manifest.json").read_text())
    for name in OUTPUT_FILES:
        assert (out / name).is_file(), name
        assert (
            manifest["outputs_sha256"][name]
            == hashlib.sha256((out / name).read_bytes()).hexdigest()
        )
    assert manifest["software_version"] == __version__
    assert manifest["effective_configuration"]["study"]["name"].startswith("AMR-Surveil synthetic")
    assert {"python", "pandas"} <= set(manifest["environment"])
    for key in ("samples", "ast", "screening", "breakpoints", "gene_drug_map"):
        assert len(manifest["inputs"][key]["sha256"]) == 64


def test_planted_association_is_recovered(demo):
    """Ground truth: P(R|gene)=0.92 vs 0.06. A working pipeline must show a large PPV/NPV gap."""
    conc = pd.read_csv(demo["out"] / "genotype_phenotype_concordance.csv")
    assert set(conc["antimicrobial"]) == set(GENE_FOR_DRUG)
    assert PLANTED["p_resistant_if_gene"] > 0.9
    pooled = conc.sum(numeric_only=True)
    ppv = pooled.n_gene_pos_R / (pooled.n_gene_pos_R + pooled.n_gene_pos_notR)
    fpr = pooled.n_gene_neg_R / (pooled.n_gene_neg_R + pooled.n_gene_neg_notR)
    assert ppv > 0.8 and fpr < 0.15


def test_planted_method_change_is_flagged_and_stratified(demo):
    flags = pd.read_csv(demo["out"] / "comparability_flags.csv")
    hit = flags[
        (flags.flag == "ast_method_changes_across_years") & (flags.antimicrobial == "Cefotaxime")
    ]
    assert len(hit) >= 1 and (hit.organism == "Escherichia coli").all()
    ind = pd.read_csv(demo["out"] / "resistance_by_group.csv")
    ecoli_cef = ind[(ind.organism == "Escherichia coli") & (ind.antimicrobial == "Cefotaxime")]
    assert {"broth_microdilution", "disk_diffusion"} <= set(
        ecoli_cef.method
    )  # never pooled across methods


def test_recorded_categories_all_rederive_from_synthetic_breakpoints(demo):
    recheck = pd.read_csv(demo["out"] / "ast_category_recheck.csv")
    assert (
        set(recheck["status"]) == {"concordant"}
        and len(recheck) == demo["summary"]["n_ast_records"]
    )


def test_heuristic_calls_excluded_by_default_and_logged(demo):
    s = demo["summary"]
    assert s["n_amr_calls"] > s["n_amr_calls_used"]
    log = pd.read_csv(demo["out"] / "exclusion_log.csv")
    assert (log.reason == "excluded_evidence_tier").sum() == s["n_amr_calls"] - s[
        "n_amr_calls_used"
    ]


def test_evidence_tier_config_changes_the_result(tmp_path):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=60)
    cfg = {**DEMO_CONFIG, "analysis": {**DEMO_CONFIG["analysis"], "excluded_evidence_tiers": []}}
    (tmp_path / "c.yaml").write_text(yaml.safe_dump(cfg))
    s = run_pipeline(
        **pipeline_kwargs(paths, config_path=tmp_path / "c.yaml", output_dir=tmp_path / "o")
    )
    assert s["n_amr_calls"] == s["n_amr_calls_used"]


def test_predictions_stay_separate_and_counts_only_by_default(demo):
    pv = pd.read_csv(demo["out"] / "prediction_vs_observed.csv")
    assert len(pv) == 3 and pv["concordance_pct"].isna().all()
    assert not (demo["out"] / "ast_with_predictions.csv").exists()


def test_attesting_independence_unlocks_metrics(tmp_path):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=120)
    cfg = {
        **DEMO_CONFIG,
        "study": {**DEMO_CONFIG["study"], "predictions_independence_attested": True},
    }
    (tmp_path / "c.yaml").write_text(yaml.safe_dump(cfg))
    run_pipeline(
        **pipeline_kwargs(paths, config_path=tmp_path / "c.yaml", output_dir=tmp_path / "o")
    )
    pv = pd.read_csv(tmp_path / "o" / "prediction_vs_observed.csv")
    assert pv["concordance_pct"].notna().all()


def test_no_go_when_study_unattested_and_no_concordance_computed(tmp_path):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=40)
    s = run_pipeline(
        **pipeline_kwargs(paths, output_dir=tmp_path / "o")
    )  # default config: not attested
    assert s["genotype_phenotype_decision"] == "NO-GO" and s["n_concordance_strata"] == 0
    ready = json.loads((tmp_path / "o" / "readiness_assessment.json").read_text())
    assert any(c["id"] == "sample_selection_documented" and not c["met"] for c in ready["criteria"])
    assert "not** calculated" in (tmp_path / "o" / "report.md").read_text()


def test_deduplication_rule_reduces_isolates_and_is_logged(tmp_path):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=120)
    base = run_pipeline(**pipeline_kwargs(paths, output_dir=tmp_path / "a"))
    cfg = {"analysis": {"deduplication_rule": "first_per_patient_species_year"}}
    (tmp_path / "c.yaml").write_text(yaml.safe_dump(cfg))
    run_pipeline(
        **pipeline_kwargs(paths, config_path=tmp_path / "c.yaml", output_dir=tmp_path / "b")
    )
    a = pd.read_csv(tmp_path / "a" / "resistance_by_group.csv").n_tested.sum()
    b = pd.read_csv(tmp_path / "b" / "resistance_by_group.csv").n_tested.sum()
    assert b < a
    log = pd.read_csv(tmp_path / "b" / "deduplication_log.csv")
    assert (log.decision == "removed").sum() > 0
    assert base["n_samples"] == 120


def test_run_is_deterministic_apart_from_timestamp(tmp_path):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=50)
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump(DEMO_CONFIG))
    for name in ("a", "b"):
        run_pipeline(**pipeline_kwargs(paths, config_path=cfg, output_dir=tmp_path / name))
    for f in OUTPUT_FILES:
        if f in {"analysis_summary.json", "report.md"}:
            continue  # carry the run timestamp
        assert (tmp_path / "a" / f).read_bytes() == (tmp_path / "b" / f).read_bytes(), f


def test_same_seed_same_data_different_seed_different_data(tmp_path):
    a = generate_synthetic_dataset(tmp_path / "a", n_samples=30, seed=1)
    b = generate_synthetic_dataset(tmp_path / "b", n_samples=30, seed=1)
    c = generate_synthetic_dataset(tmp_path / "c", n_samples=30, seed=2)
    assert a["ast"].read_bytes() == b["ast"].read_bytes() != c["ast"].read_bytes()
    with pytest.raises(ValueError):
        generate_synthetic_dataset(tmp_path / "d", n_samples=3)


def test_injected_defects_stop_strict_run_and_show_in_exploratory(tmp_path):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=60, inject_defects=True)
    with pytest.raises(QualityGateError) as err:
        run_pipeline(**pipeline_kwargs(paths, output_dir=tmp_path / "s"))
    for code in (
        "duplicate_sample_id",
        "future_collection_date",
        "invalid_ast_category",
        "duplicate_aggregate_key",
    ):
        assert code in str(err.value)
    s = run_pipeline(
        **pipeline_kwargs(
            paths, output_dir=tmp_path / "x", config_overrides={"quality": {"strict": False}}
        )
    )
    assert s["n_qa_errors"] == 4
    report = (tmp_path / "x" / "report.md").read_text()
    assert "unresolved ERROR" in report
    # the duplicated manifest ID and the orphan AST row never reach the indicators
    log = pd.read_csv(tmp_path / "x" / "exclusion_log.csv")
    assert {"ambiguous_sample_id", "unmatched_sample_id", "missing_or_invalid_result"} <= set(
        log.reason
    )
    harm = pd.read_csv(tmp_path / "x" / "label_harmonisation.csv")
    assert (harm.canonical == "Escherichia coli").any()


def test_overrides_validate_and_do_not_leak(tmp_path):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=20)
    with pytest.raises(ConfigurationError):
        run_pipeline(
            **pipeline_kwargs(
                paths,
                output_dir=tmp_path / "o",
                config_overrides={"analysis": {"confidence_level": 7}},
            )
        )
    with pytest.raises(ConfigurationError):
        run_pipeline(
            **pipeline_kwargs(
                paths, output_dir=tmp_path / "o", config_overrides={"nonsense": {"a": 1}}
            )
        )


def test_minimal_run_with_only_samples(tmp_path):
    p = tmp_path / "s.csv"
    pd.DataFrame({"sample_id": ["A", "B"]}).to_csv(p, index=False)
    s = run_pipeline(samples_path=p, output_dir=tmp_path / "o")
    assert (
        s["n_samples"] == 2
        and s["n_ast_records"] == 0
        and s["genotype_phenotype_decision"] == "NO-GO"
    )
    assert (tmp_path / "o" / "report.md").read_text().startswith("# AMR-Surveil")


def test_inventory_and_missingness_outputs(demo):
    inv = pd.read_csv(demo["out"] / "data_inventory.csv")
    assert len(inv) == 120 and inv["amr_screened"].all() and inv["has_genotype_and_ast"].all()
    assert (inv["n_ast_records"] == 3).all()
    miss = pd.read_csv(demo["out"] / "ast_missingness.csv")
    assert miss["n_missing_or_invalid_result"].sum() == 0 and miss["n_records"].sum() == 360


def test_report_has_all_sections_and_no_unrendered_placeholders(demo):
    text = (demo["out"] / "report.md").read_text()
    for heading in (
        "## Study definition",
        "## Dataset inventory",
        "## Data-quality review",
        "## Isolate-level AST resistance indicators",
        "## Comparability flags",
        "## AST category re-check",
        "## Genotype–phenotype analysis",
        "## Model predictions versus observed AST",
        "## Missingness",
        "## Interpretation checklist",
    ):
        assert heading in text, heading
    for bad in ("{", "<NA>", "nan%", "NaN", "pd.NA"):
        assert bad not in text, bad
