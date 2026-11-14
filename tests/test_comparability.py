import pandas as pd

from amr_surveil.comparability import comparability_flags


def frame(spec):
    """spec: list of (year, n, method, lab). One drug/organism/country/sector series."""
    rows, k = [], 0
    for year, n, method, lab in spec:
        for _ in range(n):
            k += 1
            rows.append(
                {
                    "sample_id": f"S{k}",
                    "country": "X",
                    "source_sector": "human",
                    "organism": "E. coli",
                    "antimicrobial": "D",
                    "year": None if year is None else str(year),
                    "method": method,
                    "breakpoint_standard": "EUCAST",
                    "breakpoint_version": "14",
                    "testing_lab": lab,
                    "specimen": "blood",
                    "result": "S",
                }
            )
    return pd.DataFrame(rows)


def flags(spec, **kw):
    out = comparability_flags(frame(spec), **kw)
    return set(out["flag"]), out


def test_clean_two_year_series_has_no_warning():
    got, _ = flags([(2022, 30, "BMD", "A"), (2023, 33, "BMD", "A")])
    assert got == set()


def test_single_year_is_info_only():
    got, out = flags([(2022, 30, "BMD", "A")])
    assert got == {"single_year_only"} and set(out["severity"]) == {"INFO"}


def test_method_lab_and_volume_changes_are_flagged():
    got, _ = flags([(2022, 20, "BMD", "A"), (2023, 60, "Disk", "B")])
    assert {
        "ast_method_changes_across_years",
        "testing_laboratory_set_changes",
        "testing_volume_jump",
    } <= got


def test_gap_year_and_small_n_and_unknown_year():
    got, _ = flags([(2020, 12, "BMD", "A"), (2022, 3, "BMD", "A"), (None, 2, "BMD", "A")])
    assert {"missing_years", "small_denominator_years", "unknown_year_records"} <= got


def test_breakpoint_version_change_is_flagged():
    f = frame([(2022, 15, "BMD", "A"), (2023, 15, "BMD", "A")])
    f.loc[f.year == "2023", "breakpoint_version"] = "15"
    assert "breakpoint_version_changes_across_years" in set(comparability_flags(f)["flag"])


def test_empty_frame():
    assert comparability_flags(frame([])).empty
