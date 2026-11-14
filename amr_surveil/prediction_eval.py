"""Descriptive comparison of model predictions with observed AST (exact sample_id + drug).

Boundaries
----------
* Predictions stay a separate stream; nothing is written back into the AST table.
* EpiResNet predicts per protein-antibiotic pair. Where several predictions exist for one
  isolate/drug, no isolate-level label is invented: those pairs are excluded and logged.
* AMR-Surveil cannot verify that the observed AST isolates were absent from the model's training
  data. Accuracy-type metrics are therefore only produced when the analyst attests independence
  (``study.predictions_independence_attested``) *and* the rows come from a held-out split;
  otherwise only counts are reported.
"""

from __future__ import annotations

import pandas as pd

from amr_surveil.epidemiology import wilson_interval
from amr_surveil.exclusions import ExclusionLog
from amr_surveil.manifest import unique_sample_ids
from amr_surveil.normalize import label_key

COLUMNS = [
    "model_name",
    "model_version",
    "checkpoint_id",
    "antimicrobial",
    "evaluation_split",
    "n_pairs",
    "n_pred_R_obs_R",
    "n_pred_R_obs_notR",
    "n_pred_S_obs_R",
    "n_pred_S_obs_notR",
    "n_obs_I_in_notR",
    "concordance_pct",
    "concordance_low_pct",
    "concordance_high_pct",
    "metrics_status",
    "interpretation_boundary",
]
BOUNDARY = (
    "Descriptive agreement between a research model and recorded AST for the same isolate and drug; "
    "independence from model training data is not verified by AMR-Surveil and no clinical performance is implied."
)


def _pred_label(value: object) -> str | None:
    key = label_key(value)
    if key is None:
        return None
    if key in {"predicted_resistant", "resistant", "r"}:
        return "R"
    if key in {"predicted_susceptible", "susceptible", "s"}:
        return "S"
    return None


def prediction_vs_observed(
    samples: pd.DataFrame,
    predictions: pd.DataFrame,
    ast: pd.DataFrame,
    *,
    independence_attested: bool = False,
    min_pairs: int = 10,
    confidence_level: float = 0.95,
    exclusions: ExclusionLog | None = None,
) -> pd.DataFrame:
    if predictions.empty or ast.empty or "predicted_class" not in predictions.columns:
        return pd.DataFrame(columns=COLUMNS)
    ids = unique_sample_ids(samples)
    pred = predictions.copy()
    pred["label"] = pred["predicted_class"].map(_pred_label)
    sid_text = pred["sample_id"].astype("string")
    if exclusions is not None:
        for rid in sid_text[~sid_text.isin(ids)].dropna().unique():
            exclusions.add(
                "prediction_eval",
                "predictions",
                rid,
                "unmatched_sample_id",
                "prediction sample_id absent from (or ambiguous in) the manifest",
            )
        unreadable = pred["label"].isna() & sid_text.isin(ids)
        for sid, cls in pred.loc[unreadable, ["sample_id", "predicted_class"]].itertuples(
            index=False
        ):
            exclusions.add(
                "prediction_eval",
                "predictions",
                sid,
                "unrecognised_predicted_class",
                f"predicted_class {cls!r} is not predicted_resistant/predicted_susceptible/R/S",
            )
    pred = pred[pred["label"].notna() & sid_text.isin(ids)].copy()
    for col in ("model_name", "model_version", "checkpoint_id", "evaluation_split"):
        if col not in pred.columns:
            pred[col] = pd.NA
    pred["drug_key"] = pred["antimicrobial"].map(label_key)
    pred["sample_id"] = pred["sample_id"].astype(str)
    # A model's train and test predictions for the same isolate/drug are different rows, not duplicates.
    mkeys = [
        "sample_id",
        "drug_key",
        "model_name",
        "model_version",
        "checkpoint_id",
        "evaluation_split",
    ]
    n_per = pred.groupby(mkeys, dropna=False)["label"].transform("size")
    if exclusions is not None:
        for sid, drug in (
            pred.loc[n_per > 1, ["sample_id", "antimicrobial"]]
            .drop_duplicates()
            .itertuples(index=False)
        ):
            exclusions.add(
                "prediction_eval",
                "predictions",
                f"{sid}|{drug}",
                "multiple_predictions_per_isolate_drug",
                "several protein-level predictions; no isolate-level label is invented",
            )
    pred = pred[n_per == 1]

    obs = ast.copy()
    obs["result"] = obs["result"].astype("string").str.strip().str.upper()
    obs = obs[obs["result"].isin(["S", "I", "R"]) & obs["sample_id"].astype(str).isin(ids)].copy()
    obs["sample_id"] = obs["sample_id"].astype(str)
    obs["drug_key"] = obs["antimicrobial"].map(label_key)
    per = obs.groupby(["sample_id", "drug_key"])["result"].agg(
        lambda s: s.iloc[0] if s.nunique() == 1 else None
    )
    if exclusions is not None:
        for sid, drug in per[per.isna()].index:
            exclusions.add(
                "prediction_eval",
                "ast",
                f"{sid}|{drug}",
                "conflicting_observed_ast_for_prediction",
                "conflicting S/I/R for the same isolate/drug; no observed label chosen",
            )
    obs = per.dropna().rename("observed").reset_index()
    joined = pred.merge(obs, on=["sample_id", "drug_key"], how="inner")
    if joined.empty:
        return pd.DataFrame(columns=COLUMNS)
    joined["obs_R"] = joined["observed"] == "R"
    joined["pred_R"] = joined["label"] == "R"
    keys = ["model_name", "model_version", "checkpoint_id", "antimicrobial", "evaluation_split"]
    rows = []
    for gkey, g in joined.groupby(keys, dropna=False, sort=True):
        a = int((g.pred_R & g.obs_R).sum())
        b = int((g.pred_R & ~g.obs_R).sum())
        c = int((~g.pred_R & g.obs_R).sum())
        d = int((~g.pred_R & ~g.obs_R).sum())
        n = a + b + c + d
        held_out = str(g["evaluation_split"].iloc[0]).strip().lower() in {
            "test",
            "held_out",
            "heldout",
        }
        if not independence_attested or not held_out:
            status = (
                "counts_only: independence from training data not attested or split not held-out"
            )
            pct = low = high = float("nan")
        elif n < min_pairs:
            status = f"counts_only: fewer than {min_pairs} pairs"
            pct = low = high = float("nan")
        else:
            status = "descriptive_metrics_attested_independent"
            lo, hi = wilson_interval(a + d, n, confidence_level)
            pct, low, high = round(100 * (a + d) / n, 4), round(lo, 4), round(hi, 4)
        row = dict(zip(keys, gkey, strict=True))
        row.update(
            n_pairs=n,
            n_pred_R_obs_R=a,
            n_pred_R_obs_notR=b,
            n_pred_S_obs_R=c,
            n_pred_S_obs_notR=d,
            n_obs_I_in_notR=int((g.observed == "I").sum()),
            concordance_pct=pct,
            concordance_low_pct=low,
            concordance_high_pct=high,
            metrics_status=status,
            interpretation_boundary=BOUNDARY,
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=COLUMNS)
