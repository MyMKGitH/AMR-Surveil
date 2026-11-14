"""Seeded mutation fuzzing: arbitrary damage to inputs may be *rejected*, but must never crash.

Every outcome has to be either a completed run (with a manifest) or an ``AMRSurveilError`` with a
readable message. Any other exception type is a bug.
"""

import random

import pandas as pd
import pytest

from amr_surveil.exceptions import AMRSurveilError
from amr_surveil.pipeline import run_pipeline
from amr_surveil.synthetic import generate_synthetic_dataset
from tests.conftest import pipeline_kwargs

HOSTILE = [
    "",
    " ",
    "NA",
    "nan",
    "None",
    "-1",
    "0",
    "1e999",
    "-inf",
    "Ünïcode✓",
    "=1+1",
    "x" * 3000,
    "R;S",
    "2024-02-30",
    "0000-00-00",
    "1999-13-01",
    "<=",
    ">",
    "<0.001",
    "99999999999999999999",
    "0.5.1",
    "\t",
    '"',
    ",",
    "r",
    "I",
]
DEDUP = ["none", "first_per_patient_species_year", "first_per_patient_species_specimen_year"]


def mutate(frame: pd.DataFrame, rng: random.Random) -> pd.DataFrame:
    frame = frame.copy()
    op = rng.choice(
        ["cell", "cell", "cell", "dropcol", "dupe", "empty", "extra", "shuffle", "blankcol"]
    )
    if op == "cell" and len(frame):
        for _ in range(rng.randint(1, 8)):
            frame.iat[rng.randrange(len(frame)), rng.randrange(len(frame.columns))] = rng.choice(
                HOSTILE
            )
    elif op == "dropcol":
        frame = frame.drop(columns=[rng.choice(list(frame.columns))])
    elif op == "dupe" and len(frame):
        frame = pd.concat(
            [frame, frame.sample(min(len(frame), 4), random_state=rng.randint(0, 99))]
        )
    elif op == "empty":
        frame = frame.iloc[0:0]
    elif op == "extra":
        frame["unexpected_column"] = "zz"
    elif op == "shuffle":
        frame = frame.sample(frac=1, random_state=rng.randint(0, 99))
    elif op == "blankcol":
        frame[rng.choice(list(frame.columns))] = ""
    return frame


@pytest.fixture(scope="module")
def base_tables(tmp_path_factory):
    paths = generate_synthetic_dataset(tmp_path_factory.mktemp("fuzzbase"), n_samples=40, seed=1)
    return {k: pd.read_csv(v, dtype=str, keep_default_na=False) for k, v in paths.items()}


@pytest.mark.parametrize("case", range(60))
def test_damaged_inputs_never_crash(tmp_path, base_tables, case):
    rng = random.Random(50_000 + case)
    touched = set(rng.sample(sorted(base_tables), rng.randint(2, 8)))
    paths = {}
    for key, table in base_tables.items():
        frame = table
        if key in touched:
            for _ in range(rng.randint(1, 3)):
                frame = mutate(frame, rng)
        paths[key] = tmp_path / f"{key}.csv"
        frame.to_csv(paths[key], index=False)
    overrides = {
        "quality": {"strict": rng.random() < 0.5},
        "analysis": {"deduplication_rule": rng.choice(DEDUP)},
    }
    try:
        run_pipeline(
            **pipeline_kwargs(paths, output_dir=tmp_path / "out", config_overrides=overrides)
        )
    except AMRSurveilError:
        return
    assert (tmp_path / "out" / "run_manifest.json").is_file()


def test_records_without_sample_id_are_logged_not_fatal(tmp_path):
    """Regression: a plasmid row with no sample_id crashed sorting under pandas 3 (found by fuzzing)."""
    (tmp_path / "s.csv").write_text("sample_id,organism\nA,Escherichia coli\n", encoding="utf-8")
    (tmp_path / "p.csv").write_text("sample_id,replicon\n,IncF\nA,IncX\n,IncY\n", encoding="utf-8")
    out = tmp_path / "o"
    run_pipeline(
        samples_path=tmp_path / "s.csv",
        plasmids_path=tmp_path / "p.csv",
        output_dir=out,
        config_overrides={"quality": {"strict": False}},
    )
    log = pd.read_csv(out / "exclusion_log.csv")
    assert (log["reason"] == "missing_sample_id").sum() == 2
