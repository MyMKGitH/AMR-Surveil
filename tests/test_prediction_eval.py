import pandas as pd
import pytest

from amr_surveil.exclusions import ExclusionLog
from amr_surveil.prediction_eval import prediction_vs_observed


def data(n=12, split="test"):
    samples = pd.DataFrame({"sample_id": [f"S{i}" for i in range(n)]})
    # observed: first half R. predicted: matches except the last of each half
    obs = ["R"] * (n // 2) + ["S"] * (n - n // 2)
    pred = (
        ["predicted_resistant"] * (n // 2 - 1)
        + ["predicted_susceptible"]
        + ["predicted_susceptible"] * (n - n // 2 - 1)
        + ["predicted_resistant"]
    )
    ast = pd.DataFrame(
        {"sample_id": samples.sample_id, "antimicrobial": "Cefotaxime", "result": obs}
    )
    preds = pd.DataFrame(
        {
            "sample_id": samples.sample_id,
            "antimicrobial": "cefotaxime",
            "predicted_class": pred,
            "model_name": "M",
            "model_version": "1",
            "checkpoint_id": "c",
            "evaluation_split": split,
        }
    )
    return samples, preds, ast


def test_counts_only_unless_attested_and_held_out():
    s, p, a = data()
    r = prediction_vs_observed(s, p, a, independence_attested=False).iloc[0]
    assert (r.n_pred_R_obs_R, r.n_pred_R_obs_notR, r.n_pred_S_obs_R, r.n_pred_S_obs_notR) == (
        5,
        1,
        1,
        5,
    )
    assert pd.isna(r.concordance_pct) and r.metrics_status.startswith("counts_only")
    s, p, a = data(split="train")
    r = prediction_vs_observed(s, p, a, independence_attested=True).iloc[0]
    assert pd.isna(r.concordance_pct)  # attested but not a held-out split


def test_metrics_when_attested_and_held_out():
    s, p, a = data()
    r = prediction_vs_observed(s, p, a, independence_attested=True, min_pairs=10).iloc[0]
    assert r.concordance_pct == pytest.approx(100 * 10 / 12, abs=1e-3)
    assert r.concordance_low_pct < r.concordance_pct < r.concordance_high_pct
    assert r.metrics_status == "descriptive_metrics_attested_independent"
    small = prediction_vs_observed(s, p, a, independence_attested=True, min_pairs=50).iloc[0]
    assert pd.isna(small.concordance_pct)


def test_protein_level_duplicates_get_no_invented_isolate_label():
    s, p, a = data()
    dup = p.iloc[[0]].assign(predicted_class="predicted_susceptible")
    log = ExclusionLog()
    out = prediction_vs_observed(s, pd.concat([p, dup]), a, exclusions=log)
    assert out.iloc[0].n_pairs == 11
    assert "multiple_predictions_per_isolate_drug" in set(log.to_frame()["reason"])


def test_unobserved_or_conflicting_ast_is_not_joined():
    s, p, a = data()
    conflict = pd.concat([a, a.iloc[[0]].assign(result="S")])
    assert prediction_vs_observed(s, p, conflict).iloc[0].n_pairs == 11
    assert prediction_vs_observed(s, p, a.iloc[0:0]).empty
    assert prediction_vs_observed(s, p.assign(predicted_class="research_model_output"), a).empty
