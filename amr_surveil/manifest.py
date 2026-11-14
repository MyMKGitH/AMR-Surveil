"""Helpers for the sample manifest: only unambiguous sample keys may be joined."""

from __future__ import annotations

import pandas as pd


def unique_sample_ids(samples: pd.DataFrame) -> set[str]:
    """Sample IDs that occur exactly once in the manifest (safe for exact-ID joins)."""
    if "sample_id" not in samples.columns or samples.empty:
        return set()
    counts = samples["sample_id"].dropna().astype(str).value_counts()
    return set(counts[counts == 1].index)


def ambiguous_sample_ids(samples: pd.DataFrame) -> set[str]:
    """Sample IDs that occur more than once; joins on them would be ambiguous."""
    if "sample_id" not in samples.columns or samples.empty:
        return set()
    counts = samples["sample_id"].dropna().astype(str).value_counts()
    return set(counts[counts > 1].index)


def unique_manifest(samples: pd.DataFrame) -> pd.DataFrame:
    """The manifest restricted to unambiguous sample keys."""
    if "sample_id" not in samples.columns:
        return pd.DataFrame(columns=samples.columns)
    keep = unique_sample_ids(samples)
    return samples[samples["sample_id"].astype(str).isin(keep)].copy()
