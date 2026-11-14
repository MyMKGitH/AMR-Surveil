"""Hand-computed 2x2 tables: the numbers below are derived on paper, not from the code under test."""

import pandas as pd
import pytest

from amr_surveil.exclusions import ExclusionLog
from amr_surveil.genotype_phenotype import assess_readiness, genotype_phenotype_concordance

ORG = "Escherichia coli"
STUDY_OK = {"sample_selection_documented": True, "selection_description": "all isolates 2024"}


def build(n_tp=6, n_fp=1, n_fn=2, n_tn=11, *, screened_all=True):
    """TP gene+/R, FP gene+/notR, FN gene-/R, TN gene-/notR (one of the TN is I)."""
    rows, calls, ast, scr = [], [], [], []
    k = 0
    for count, gene_pos, res in (
        (n_tp, True, "R"),
        (n_fp, True, "S"),
        (n_fn, False, "R"),
        (n_tn, False, "S"),
    ):
        for _ in range(count):
            k += 1
            sid = f"I{k:03d}"
            rows.append({"sample_id": sid, "organism": ORG, "collection_date": "2024-05-01"})
            if gene_pos:
                calls.append(
                    {"sample_id": sid, "gene": "blaCTX-M-15", "evidence_tier": "curated_database"}
                )
            result = "I" if (res == "S" and not gene_pos and k == n_tp + n_fp + n_fn + 1) else res
            ast.append(
                {
                    "sample_id": sid,
                    "organism": ORG,
                    "antimicrobial": "Cefotaxime",
                    "result": result,
                    "method": "BMD",
                    "breakpoint_standard": "SYN",
                    "breakpoint_version": "1",
                }
            )
            scr.append(
                {
                    "sample_id": sid,
                    "method": "AMRFinderPlus",
                    "database": "DB",
                    "database_version": "v1",
                }
            )
    gm = pd.DataFrame(
        [
            {
                "organism": ORG,
                "gene": "blaCTX-M-15",
                "antimicrobial": "Cefotaxime",
                "mapping_source": "Fixture citation",
                "match_type": "exact",
            }
        ]
    )
    return (
        pd.DataFrame(rows),
        pd.DataFrame(calls, columns=["sample_id", "gene", "evidence_tier"]),
        pd.DataFrame(ast),
        pd.DataFrame(scr),
        gm,
    )


def test_concordance_matches_hand_computation():
    samples, calls, ast, scr, gm = build()
    summary, disc = genotype_phenotype_concordance(samples, calls, ast, scr, gm, min_isolates=10)
    r = summary.iloc[0]
    assert (r.n_gene_pos_R, r.n_gene_pos_notR, r.n_gene_neg_R, r.n_gene_neg_notR) == (6, 1, 2, 11)
    assert r.n_I_in_notR == 1  # I is not-R but reported separately
    assert r.sensitivity_pct == pytest.approx(100 * 6 / 8)
    assert r.specificity_pct == pytest.approx(100 * 11 / 12)
    assert r.ppv_pct == pytest.approx(100 * 6 / 7)
    assert r.npv_pct == pytest.approx(100 * 11 / 13)
    assert not r.metrics_suppressed
    assert r.sensitivity_low_pct < r.sensitivity_pct < r.sensitivity_high_pct
    assert set(disc["discordance_type"]) == {
        "mapped_gene_detected_but_phenotype_not_R",
        "phenotype_R_without_mapped_gene",
    }
    assert len(disc) == 3


def test_small_strata_report_counts_only():
    samples, calls, ast, scr, gm = build(2, 0, 1, 3)
    summary, _ = genotype_phenotype_concordance(samples, calls, ast, scr, gm, min_isolates=10)
    r = summary.iloc[0]
    assert r.n_isolates == 6 and r.metrics_suppressed
    assert pd.isna(r.sensitivity_pct) and pd.isna(r.ppv_low_pct)


def test_unscreened_isolates_are_never_negatives():
    samples, calls, ast, scr, gm = build()
    kept = scr.iloc[:10]  # only 10 of 20 screened
    summary, _ = genotype_phenotype_concordance(samples, calls, ast, kept, gm, min_isolates=1)
    assert summary["n_isolates"].sum() == 10


def test_conflicting_screening_records_exclude_the_isolate():
    samples, calls, ast, scr, gm = build()
    extra = scr.iloc[[0]].assign(database_version="v2")
    log = ExclusionLog()
    summary, _ = genotype_phenotype_concordance(
        samples, calls, ast, pd.concat([scr, extra]), gm, min_isolates=1, exclusions=log
    )
    assert summary["n_isolates"].sum() == 19
    assert "multiple_screening_records" in set(log.to_frame()["reason"])


def test_gene_for_other_drug_does_not_count_as_positive():
    samples, calls, ast, scr, gm = build()
    gm2 = gm.assign(gene="blaKPC-2")  # mapping no longer links the carried gene to this drug
    summary, _ = genotype_phenotype_concordance(samples, calls, ast, scr, gm2, min_isolates=1)
    assert summary.iloc[0].n_gene_pos_R + summary.iloc[0].n_gene_pos_notR == 0


def test_prefix_match_type_and_exact_default():
    samples, calls, ast, scr, gm = build()
    calls = calls.assign(gene="blaCTX-M-27")
    exact, _ = genotype_phenotype_concordance(samples, calls, ast, scr, gm, min_isolates=1)
    assert exact.iloc[0].n_gene_pos_R == 0
    prefix, _ = genotype_phenotype_concordance(
        samples, calls, ast, scr, gm.assign(gene="blaCTX-M", match_type="prefix"), min_isolates=1
    )
    assert prefix.iloc[0].n_gene_pos_R == 6


def test_organism_mismatch_between_manifest_and_ast_is_excluded():
    samples, calls, ast, scr, gm = build()
    ast.loc[0, "organism"] = "Klebsiella pneumoniae"
    log = ExclusionLog()
    summary, _ = genotype_phenotype_concordance(
        samples, calls, ast, scr, gm, min_isolates=1, exclusions=log
    )
    assert summary["n_isolates"].sum() == 19


def test_readiness_go_and_each_blocker():
    samples, calls, ast, scr, gm = build()
    assert assess_readiness(samples, calls, ast, scr, gm, STUDY_OK)["decision"] == "GO"

    def blocked(**kw):
        args = dict(
            samples=samples, amr_calls=calls, ast=ast, screening=scr, gene_map=gm, study=STUDY_OK
        )
        args.update(kw)
        res = assess_readiness(**args)
        assert res["decision"] == "NO-GO"
        return {c["id"] for c in res["criteria"] if not c["met"]}

    assert "sample_selection_documented" in blocked(study={})
    assert "cited_gene_drug_map" in blocked(gene_map=gm.assign(mapping_source=pd.NA))
    assert "cited_gene_drug_map" in blocked(gene_map=gm.iloc[0:0])
    assert "ast_method_recorded" in blocked(ast=ast.assign(method=pd.NA))
    assert "interpretive_standard_recorded" in blocked(ast=ast.assign(breakpoint_version=pd.NA))
    assert "screening_provenance_recorded" in blocked(screening=scr.assign(database_version=pd.NA))
    assert "paired_isolate_level" in blocked(screening=scr.iloc[0:0])
