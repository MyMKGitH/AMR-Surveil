import math

import pandas as pd
import pytest

from amr_surveil.normalize import harmonise_labels, label_key, parse_iso_dates, parse_mic


def test_iso_dates_precision_and_rejections():
    out = parse_iso_dates(
        pd.Series(["2024-03-05", "2024-03", "2024", "03/05/2024", "2024-02-30", "20240305", None])
    )
    assert out["year"].tolist()[:3] == [2024, 2024, 2024]
    assert out["precision"].tolist()[:3] == ["day", "month", "year"]
    assert bool(out["non_iso"].iloc[3])  # ambiguous day/month order is never guessed
    assert bool(out["impossible"].iloc[4])  # 30 Feb
    assert out["year"].iloc[3] is pd.NA or pd.isna(out["year"].iloc[3])
    assert pd.isna(out["year"].iloc[6])


@pytest.mark.parametrize(
    "raw,op,val,ok",
    [
        ("0.25", "=", 0.25, True),
        ("<=0.5", "<=", 0.5, True),
        (">32", ">", 32.0, True),
        ("≥8", ">=", 8.0, True),
        ("≤ 2", "<=", 2.0, True),
        ("1,5", None, None, False),
        ("abc", None, None, False),
        ("", None, None, False),
        (None, None, None, False),
        ("-1", None, None, False),
    ],
)
def test_parse_mic(raw, op, val, ok):
    m = parse_mic(raw)
    assert m.valid is ok
    if ok:
        assert m.operator == op and m.value == val
    else:
        assert math.isnan(m.value)


def test_harmonise_labels_folds_case_not_synonyms():
    s = pd.Series(["Escherichia coli", "escherichia  coli ", "Escherichia coli", "E. coli", None])
    out, merged = harmonise_labels(s)
    assert out.iloc[1] == "Escherichia coli"
    assert out.iloc[3] == "E. coli"  # synonyms need a curated dictionary; never guessed
    assert pd.isna(out.iloc[4])
    assert len(merged) == 1 and merged[0]["canonical"] == "Escherichia coli"


def test_label_key_handles_missing():
    assert label_key(None) is None and label_key(pd.NA) is None
    assert label_key("  A   b ") == "a b"


def test_write_csv_neutralises_spreadsheet_formulas(tmp_path):
    import pandas as pd

    from amr_surveil.io import write_csv

    frame = pd.DataFrame(
        {
            "sample_id": ["=HYPERLINK(\"http://x\")", "+cmd", "@SUM(1)", "-cmd", "-1.5", "S1", None],
            "n": [1, 2, 3, 4, 5, 6, 7],
            "v": [-1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        }
    )
    target = tmp_path / "o.csv"
    write_csv(frame, target)
    text = target.read_text(encoding="utf-8").splitlines()
    assert text[1].startswith("\"'=HYPERLINK")
    assert text[2].startswith("'+cmd") and text[3].startswith("'@SUM")
    assert text[4].startswith("'-cmd")
    assert text[5].startswith("-1.5,") and text[6].startswith("S1,")
    assert ",-1.0" in text[1]  # numeric columns untouched
