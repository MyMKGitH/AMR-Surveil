"""Second, independent implementation (stdlib only, no pandas) of the headline numbers.

It re-derives resistance proportions, Wilson intervals and genotype-phenotype 2x2 counts directly
from the raw CSV inputs and compares them with the pipeline outputs. Agreement between two
implementations written differently is evidence the counting rules are implemented as documented.
"""

import collections
import csv
import math

import pytest

from amr_surveil.pipeline import run_pipeline
from amr_surveil.synthetic import generate_synthetic_dataset
from tests.conftest import pipeline_kwargs

CONFIG = {
    "analysis": {"reference_date": "2026-11-07"},
    "quality": {"strict": True},
    "study": {
        "name": "oracle",
        "question": "q",
        "inclusion_criteria": "all",
        "exclusion_criteria": "none",
        "sample_selection_documented": True,
        "selection_description": "synthetic",
    },
}
WEAK_TIERS = ("heuristic_motif", "homology_search")


def read_rows(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def wilson(k, n, z=1.959963984540054):
    p = k / n
    denom = 1 + z * z / n
    centre = p + z * z / (2 * n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (centre - half) / denom * 100, (centre + half) / denom * 100


def eligible(ast, samples):
    for r in ast:
        s = samples.get(r["sample_id"])
        res = r["result"].strip().upper()
        if s is None or res not in ("S", "I", "R"):
            continue
        if r["organism"] and r["organism"].strip() != s["organism"]:
            continue
        yield s, r, res


@pytest.mark.parametrize(("seed", "n"), [(1, 60), (2, 120), (3, 300), (5, 40), (20261107, 120)])
def test_pipeline_matches_independent_recomputation(tmp_path, seed, n):
    paths = generate_synthetic_dataset(tmp_path / "in", n_samples=n, seed=seed)
    out = tmp_path / "out"
    run_pipeline(**pipeline_kwargs(paths, output_dir=out, config_overrides=CONFIG))
    samples = {r["sample_id"]: r for r in read_rows(paths["samples"])}
    ast = read_rows(paths["ast"])

    # ---- resistance proportions
    groups = collections.defaultdict(lambda: [0, 0, 0, 0])
    for s, r, res in eligible(ast, samples):
        key = (
            s["collection_date"][:4],
            s["country"],
            s["source_sector"],
            s["organism"],
            r["antimicrobial"],
            r["method"],
            r["breakpoint_standard"],
            r["breakpoint_version"],
        )
        g = groups[key]
        g[0] += 1
        g[1] += res == "R"
        g[2] += res == "S"
        g[3] += res == "I"
    got = read_rows(out / "resistance_by_group.csv")
    assert len(got) == len(groups)
    for r in got:
        key = tuple(
            r[c]
            for c in (
                "year",
                "country",
                "source_sector",
                "organism",
                "antimicrobial",
                "method",
                "breakpoint_standard",
                "breakpoint_version",
            )
        )
        tested, resistant, s_n, i_n = groups[key]
        assert (int(r["n_tested"]), int(r["n_resistant"]), int(r["n_S"]), int(r["n_I"])) == (
            tested,
            resistant,
            s_n,
            i_n,
        )
        low, high = wilson(resistant, tested)
        assert float(r["wilson_low_pct"]) == pytest.approx(low, abs=1e-3)
        assert float(r["wilson_high_pct"]) == pytest.approx(high, abs=1e-3)
        assert float(r["resistant_pct"]) == pytest.approx(100 * resistant / tested, abs=1e-3)

    # ---- genotype-phenotype 2x2
    screening = read_rows(paths["screening"])
    freq = collections.Counter(r["sample_id"] for r in screening)
    screened = {r["sample_id"]: r for r in screening if freq[r["sample_id"]] == 1}
    genes = collections.defaultdict(set)
    for r in read_rows(paths["amr_calls"]):
        if r["evidence_tier"] not in WEAK_TIERS:
            genes[r["sample_id"]].add(r["gene"])
    gene_map = collections.defaultdict(set)
    for r in read_rows(paths["gene_drug_map"]):
        gene_map[(r["organism"], r["antimicrobial"])].add(r["gene"])
    counts = collections.defaultdict(lambda: [0, 0, 0, 0, 0])
    for s, r, res in eligible(ast, samples):
        if r["sample_id"] not in screened:
            continue
        screen = screened[r["sample_id"]]
        positive = bool(genes[r["sample_id"]] & gene_map[(s["organism"], r["antimicrobial"])])
        key = (
            s["organism"],
            r["antimicrobial"],
            r["method"],
            r["breakpoint_standard"],
            r["breakpoint_version"],
            screen["database"],
            screen["database_version"],
        )
        c = counts[key]
        resistant = res == "R"
        c[0 if positive and resistant else 1 if positive else 2 if resistant else 3] += 1
        c[4] += res == "I"
    got = read_rows(out / "genotype_phenotype_concordance.csv")
    assert len(got) == len(counts)
    for r in got:
        key = tuple(
            r[c]
            for c in (
                "organism",
                "antimicrobial",
                "method",
                "breakpoint_standard",
                "breakpoint_version",
                "screen_database",
                "screen_database_version",
            )
        )
        tp, fp, fn, tn, n_i = counts[key]
        assert tuple(
            int(r[c])
            for c in (
                "n_gene_pos_R",
                "n_gene_pos_notR",
                "n_gene_neg_R",
                "n_gene_neg_notR",
                "n_I_in_notR",
            )
        ) == (tp, fp, fn, tn, n_i)
        if r["metrics_suppressed"].lower() != "true":
            assert float(r["sensitivity_pct"]) == pytest.approx(100 * tp / (tp + fn), abs=1e-3)
            assert float(r["specificity_pct"]) == pytest.approx(100 * tn / (tn + fp), abs=1e-3)
