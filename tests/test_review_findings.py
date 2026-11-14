"""Regression tests for defects found by an independent adversarial review (see CHANGELOG 1.0.0).

Each test reproduces the reviewer's failing input and pins the corrected behaviour.
"""

import pandas as pd
import pytest

from amr_surveil.ast_interpretation import recheck_ast_categories
from amr_surveil.epidemiology import deduplicate_isolates
from amr_surveil.exclusions import ExclusionLog
from amr_surveil.normalize import parse_mic
from amr_surveil.pipeline import run_pipeline
from amr_surveil.prediction_eval import prediction_vs_observed
from amr_surveil.synthetic import generate_synthetic_dataset
from tests.conftest import pipeline_kwargs
from tests.test_independent_oracle import CONFIG


def run_small(tmp_path, samples, ast, **extra):
    (tmp_path / "s.csv").write_text(samples, encoding="utf-8")
    (tmp_path / "a.csv").write_text(ast, encoding="utf-8")
    kwargs = {"ast_path": tmp_path / "a.csv"}
    for key, text in extra.items():
        (tmp_path / f"{key}.csv").write_text(text, encoding="utf-8")
        kwargs[f"{key}_path"] = tmp_path / f"{key}.csv"
    run_pipeline(
        samples_path=tmp_path / "s.csv",
        output_dir=tmp_path / "out",
        config_overrides={"quality": {"strict": False}, **kwargs.pop("config", {})},
        **kwargs,
    )
    return tmp_path / "out"


# 1 ---------------------------------------------------------------- year-only dates
def test_year_only_date_makes_first_isolate_order_uncertain():
    samples = pd.DataFrame(
        {
            "sample_id": ["a", "b", "c"],
            "patient_id": ["P1"] * 3,
            "organism": ["Escherichia coli"] * 3,
            "collection_date": ["2024-06-01", "2024-07-01", "2024"],
        }
    )
    _, log = deduplicate_isolates(samples, "first_per_patient_species_year")
    first = log[log.decision == "kept"].iloc[0]
    assert bool(first.date_tie) is True
    assert "uncertain" in first.reason


# 2 -------------------------------------------------- dedup must consider tested isolates
def test_untested_first_isolate_does_not_remove_tested_repeat(tmp_path):
    samples = (
        "sample_id,organism,collection_date,patient_id,country,source_sector\n"
        "a,Escherichia coli,2024-01-01,P1,NO,human\n"
        "b,Escherichia coli,2024-03-01,P1,NO,human\n"
    )
    ast = "sample_id,antimicrobial,result,method\nb,ciprofloxacin,R,BMD\n"
    out = tmp_path / "o"
    (tmp_path / "s.csv").write_text(samples)
    (tmp_path / "a.csv").write_text(ast)
    run_pipeline(
        samples_path=tmp_path / "s.csv",
        ast_path=tmp_path / "a.csv",
        output_dir=out,
        config_overrides={
            "quality": {"strict": False},
            "analysis": {"deduplication_rule": "first_per_patient_species_year"},
        },
    )
    ind = pd.read_csv(out / "resistance_by_group.csv")
    assert ind["n_tested"].sum() == 1 and ind["n_resistant"].sum() == 1


# 3 ------------------------------------------- organism spelling shared across tables
def test_organism_case_variant_across_tables_does_not_split_species(tmp_path):
    rows = "".join(f"s{i},Escherichia coli,2024-01-01,NO,human\n" for i in range(10))
    samples = "sample_id,organism,collection_date,country,source_sector\n" + rows
    ast = "sample_id,organism,antimicrobial,result,method\n" + "".join(
        f"s{i},{'escherichia coli' if i < 6 else ''},cip,{'R' if i % 2 else 'S'},BMD\n"
        for i in range(10)
    )
    out = run_small(tmp_path, samples, ast)
    ind = pd.read_csv(out / "resistance_by_group.csv")
    assert list(ind["organism"]) == ["Escherichia coli"] and int(ind["n_tested"].iloc[0]) == 10


# 4 ------------------------------------------------------------------ MIC units
BP = pd.DataFrame(
    {
        "organism": ["E. coli"],
        "antimicrobial": ["cip"],
        "standard": ["EUCAST"],
        "version": ["14"],
        "mic_s_le": [0.25],
        "mic_r_gt": [0.5],
    }
)


def one_record(**fields):
    base = {
        "sample_id": "a",
        "organism": "E. coli",
        "antimicrobial": "cip",
        "result": "R",
        "breakpoint_standard": "EUCAST",
        "breakpoint_version": "14",
    }
    return pd.DataFrame([{**base, **fields}])


@pytest.mark.parametrize(
    ("fields", "status", "derived"),
    [
        ({"mic": "1"}, "concordant", "R"),  # blank unit = mg/L
        ({"mic": "1", "mic_unit": "mg/L"}, "concordant", "R"),
        ({"mic": "0.1", "mic_unit": "µg/mL"}, "discordant", "S"),  # equivalent unit
        ({"mic": "200", "mic_unit": "ug/L"}, "unit_not_mg_per_l", None),  # 0.2 mg/L, not 200
    ],
)
def test_mic_unit_is_respected(fields, status, derived):
    out = recheck_ast_categories(one_record(**fields), BP).iloc[0]
    assert out.status == status
    assert (out.derived_result if pd.notna(out.derived_result) else None) == derived


# 5 ------------------------------------------- unknown standard/version never matches
def test_missing_standard_or_version_never_matches_a_breakpoint_row():
    bp = BP.assign(standard=None, version=None)
    rec = one_record(mic="0.1", breakpoint_standard=None, breakpoint_version=None)
    assert recheck_ast_categories(rec, bp).iloc[0].status == "no_breakpoint"


# 6 --------------------------------- organism missing in AST is filled from the manifest
def test_recheck_uses_manifest_organism_when_ast_organism_is_blank(tmp_path):
    samples = (
        "sample_id,organism,collection_date,country,source_sector\nA,E. coli,2024-01-01,NO,human\n"
    )
    ast = (
        "sample_id,antimicrobial,result,method,mic,breakpoint_standard,breakpoint_version\n"
        "A,cip,R,BMD,4,EUCAST,14\n"
    )
    bp = "organism,antimicrobial,standard,version,mic_s_le,mic_r_gt\nE. coli,cip,EUCAST,14,0.25,0.5\n"
    out = run_small(tmp_path, samples, ast, breakpoints=bp)
    recheck = pd.read_csv(out / "ast_category_recheck.csv")
    assert recheck.iloc[0].status == "concordant" and recheck.iloc[0].derived_result == "R"


# 7 -------------------------------------- weak / partial hits are not genotype-negative
def _gp_inputs(tmp_path, mutate_calls):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=80, seed=7)
    calls = pd.read_csv(paths["amr_calls"], dtype=str, keep_default_na=False)
    mutate_calls(calls).to_csv(paths["amr_calls"], index=False)
    run_pipeline(**pipeline_kwargs(paths, output_dir=tmp_path / "out", config_overrides=CONFIG))
    return pd.read_csv(tmp_path / "out" / "genotype_phenotype_concordance.csv")


def _ctx_carriers(tmp_path_factory_dir):
    paths = generate_synthetic_dataset(tmp_path_factory_dir / "probe", n_samples=80, seed=7)
    calls = pd.read_csv(paths["amr_calls"], dtype=str, keep_default_na=False)
    return calls[
        (calls.gene == "blaCTX-M-15") & (calls.organism == "Escherichia coli")
    ].sample_id.tolist()


@pytest.mark.parametrize(
    ("column", "value"), [("evidence_tier", "homology_search"), ("hit_method", "PARTIALX")]
)
def test_low_evidence_hit_is_excluded_not_counted_negative(tmp_path, column, value):
    carriers = _ctx_carriers(tmp_path)
    assert len(carriers) >= 2
    chosen = carriers[:2]

    baseline = _gp_inputs(tmp_path / "base", lambda c: c)
    changed = _gp_inputs(
        tmp_path / "chg",
        lambda c: (
            c.assign(**{column: c[column].where(~c.sample_id.isin(chosen), value)})
            if column in c
            else c
        ),
    )
    pick = lambda df: df[(df.organism == "Escherichia coli") & (df.antimicrobial == "Cefotaxime")]  # noqa: E731
    base_rows, new_rows = pick(baseline), pick(changed)
    cols = ["n_gene_pos_R", "n_gene_pos_notR", "n_gene_neg_R", "n_gene_neg_notR"]
    assert int(new_rows["n_excluded_uncertain_genotype"].sum()) >= 1
    # excluded isolates vanish from the 2x2 table; they do NOT show up as extra negatives
    assert new_rows[cols].to_numpy().sum() < base_rows[cols].to_numpy().sum()
    assert new_rows["n_gene_neg_R"].sum() <= base_rows["n_gene_neg_R"].sum()


# 8 ---------------------------------------------- nothing dropped silently from predictions
def test_prediction_drops_are_logged():
    samples = pd.DataFrame({"sample_id": ["s1", "s2"], "organism": ["E"] * 2})
    pred = pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "ghost"],
            "antimicrobial": ["cip"] * 3,
            "predicted_class": ["1", "predicted_resistant", "predicted_resistant"],
            "evaluation_split": ["test"] * 3,
        }
    )
    ast = pd.DataFrame(
        {
            "sample_id": ["s2", "s2"],
            "antimicrobial": ["cip", "cip"],
            "result": ["R", "S"],  # conflicting observed results
        }
    )
    log = ExclusionLog()
    prediction_vs_observed(samples, pred, ast, exclusions=log)
    reasons = set(log.to_frame()["reason"])
    assert {
        "unrecognised_predicted_class",
        "unmatched_sample_id",
        "conflicting_observed_ast_for_prediction",
    } <= reasons


def test_train_and_test_predictions_for_one_isolate_are_not_duplicates():
    samples = pd.DataFrame({"sample_id": ["s1"], "organism": ["E"]})
    pred = pd.DataFrame(
        {
            "sample_id": ["s1", "s1"],
            "antimicrobial": ["cip", "cip"],
            "predicted_class": ["predicted_resistant"] * 2,
            "evaluation_split": ["train", "test"],
        }
    )
    ast = pd.DataFrame({"sample_id": ["s1"], "antimicrobial": ["cip"], "result": ["R"]})
    out = prediction_vs_observed(samples, pred, ast)
    assert sorted(out["evaluation_split"]) == ["test", "train"]


# 9 ------------------------------------------------------------ blank antimicrobial
def test_blank_antimicrobial_is_excluded_not_pooled_as_unspecified(tmp_path):
    samples = "sample_id,organism,collection_date,country,source_sector\nA,E,2024-01-01,NO,human\n"
    ast = "sample_id,antimicrobial,result,method\nA,,R,BMD\nA,cip,S,BMD\n"
    out = run_small(tmp_path, samples, ast)
    ind = pd.read_csv(out / "resistance_by_group.csv")
    assert list(ind["antimicrobial"]) == ["cip"]
    assert "missing_antimicrobial" in set(pd.read_csv(out / "exclusion_log.csv")["reason"])


# minor ------------------------------------------------------------------------
@pytest.mark.parametrize("text", ["1.", "16.", ".5", "0.5"])
def test_mic_with_trailing_or_leading_dot_parses(text):
    assert parse_mic(text).valid


def test_concordance_honours_deduplication(tmp_path):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=120, seed=3)
    cfg = {**CONFIG, "analysis": {**CONFIG["analysis"], "deduplication_rule": "none"}}
    run_pipeline(**pipeline_kwargs(paths, output_dir=tmp_path / "none", config_overrides=cfg))
    cfg["analysis"] = {**cfg["analysis"], "deduplication_rule": "first_per_patient_species_year"}
    run_pipeline(**pipeline_kwargs(paths, output_dir=tmp_path / "dedup", config_overrides=cfg))
    n_none = pd.read_csv(tmp_path / "none" / "genotype_phenotype_concordance.csv")[
        "n_isolates"
    ].sum()
    n_dedup = pd.read_csv(tmp_path / "dedup" / "genotype_phenotype_concordance.csv")[
        "n_isolates"
    ].sum()
    assert n_dedup < n_none
