from __future__ import annotations

from amr_surveil import doctor
from amr_surveil.cli import main


def test_doctor_reports_missing_optional_module_as_warning() -> None:
    check = doctor.check_module(
        "definitely_not_installed_xyz", required=False, purpose="x", fix="pip install it"
    )
    assert check.status == doctor.WARN and check.fix


def test_doctor_reports_missing_required_module_as_failure() -> None:
    check = doctor.check_module("definitely_not_installed_xyz", required=True, purpose="x", fix="f")
    assert check.status == doctor.FAIL


def test_doctor_catches_broken_binary_import(monkeypatch) -> None:
    def boom(name: str):
        raise ImportError("DLL load failed")

    monkeypatch.setattr(doctor.importlib, "import_module", boom)
    version, error = doctor._import_version("pyarrow")
    assert version is None and "DLL load failed" in error


def test_doctor_smoke_run_and_cli(capsys) -> None:
    assert doctor.check_smoke().status == doctor.OK
    code = main(["doctor", "--no-smoke"])
    out = capsys.readouterr().out
    assert code in (0, 1) and "python" in out


def test_exit_code_and_render() -> None:
    checks = [doctor.Check("a", doctor.OK, "1"), doctor.Check("b", doctor.FAIL, "bad", "do x")]
    assert doctor.exit_code(checks) == 1
    assert "fix: do x" in doctor.render(checks)
