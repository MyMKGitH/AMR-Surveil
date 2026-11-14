import pandas as pd
import pytest

from amr_surveil.ast_interpretation import classify_mic, classify_zone, recheck_ast_categories

NAN = float("nan")


@pytest.mark.parametrize(
    "op,val,expected",
    [
        ("=", 1, "S"),
        ("=", 1.5, "I"),
        ("=", 2, "I"),
        ("=", 4, "R"),
        ("<=", 1, "S"),
        ("<=", 2, None),
        ("<", 0.25, "S"),
        (">", 4, "R"),
        (">", 2, "R"),
        (">", 1.5, None),
        ("<", 2, None),
        (">=", 4, "R"),
        (">=", 2, None),
        ("=", NAN, None),
    ],
)
def test_classify_mic_with_censoring(op, val, expected):
    # synthetic breakpoint: S <= 1, R > 2
    assert classify_mic(op, val, 1.0, 2.0) == expected


def test_classify_zone():
    assert classify_zone(20, 20, 17) == "S"
    assert classify_zone(18, 20, 17) == "I"
    assert classify_zone(16, 20, 17) == "R"
    assert classify_zone(NAN, 20, 17) is None


def make(rows):
    base = {
        "sample_id": "A",
        "organism": "E. coli",
        "antimicrobial": "Cefotaxime",
        "result": "S",
        "mic": "0.5",
        "zone_mm": None,
        "breakpoint_standard": "SYN",
        "breakpoint_version": "1",
    }
    return pd.DataFrame([{**base, **r} for r in rows])


BP = pd.DataFrame(
    [
        {
            "organism": "e. coli",
            "antimicrobial": "cefotaxime",
            "standard": "syn",
            "version": "1",
            "mic_s_le": "1",
            "mic_r_gt": "2",
            "zone_s_ge": "20",
            "zone_r_lt": "17",
        }
    ]
)


def test_recheck_statuses():
    ast = make(
        [
            {"sample_id": "ok"},
            {"sample_id": "bad", "result": "R"},
            {"sample_id": "cens", "mic": ">0.5", "result": "S"},  # indeterminate
            {"sample_id": "zone", "mic": None, "zone_mm": "15", "result": "R"},  # zone path
            {"sample_id": "nomeas", "mic": None},
            {"sample_id": "nobp", "breakpoint_version": "2"},  # version not in table
        ]
    )
    out = recheck_ast_categories(ast, BP).set_index("sample_id")["status"].to_dict()
    assert out == {
        "ok": "concordant",
        "bad": "discordant",
        "cens": "indeterminate",
        "zone": "concordant",
        "nomeas": "no_measurement",
        "nobp": "no_breakpoint",
    }


def test_ambiguous_breakpoint_rows_are_never_guessed():
    twice = pd.concat([BP, BP.assign(mic_s_le="0.1")])
    out = recheck_ast_categories(make([{}]), twice)
    assert out.iloc[0]["status"] == "no_breakpoint"


def test_no_breakpoints_supplied_returns_empty():
    assert recheck_ast_categories(make([{}]), BP.iloc[0:0]).empty


def test_blank_string_mic_unit_is_assumed_mg_per_l_and_wilson_boundaries_are_exact():
    from amr_surveil.ast_interpretation import _mic_unit_ok
    from amr_surveil.epidemiology import wilson_interval

    assert _mic_unit_ok("") and _mic_unit_ok("  ") and _mic_unit_ok(None)
    assert not _mic_unit_ok("nM")
    assert wilson_interval(0, 10)[0] == 0.0
    assert wilson_interval(10, 10)[1] == 100.0
    low, high = wilson_interval(0, 10)
    assert abs(high - 27.7533) < 1e-3  # published Wilson value for 0/10
