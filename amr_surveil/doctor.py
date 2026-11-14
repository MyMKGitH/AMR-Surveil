"""Environment self-check (``amr-surveil doctor``).

Most first-run failures are environment problems, not analysis problems (a broken binary wheel, an
interpreter newer than the dependency wheels, a missing optional UI package). This module reports
them in plain words with a suggested fix, and finishes with an end-to-end smoke run on a tiny
synthetic dataset so "OK" means the pipeline actually executes here.
"""

from __future__ import annotations

import importlib
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

OK, WARN, FAIL = "OK", "WARN", "FAIL"


@dataclass(frozen=True)
class Check:
    name: str
    status: str
    detail: str
    fix: str = ""


def _import_version(module: str) -> tuple[str | None, str]:
    """Return (version, error). A module that imports but has no ``__version__`` reports ``unknown``."""
    try:
        mod = importlib.import_module(module)
    except Exception as exc:  # broken binary wheels raise ImportError subclasses or worse
        return None, f"{type(exc).__name__}: {exc}"
    return str(getattr(mod, "__version__", "unknown")), ""


def check_python() -> Check:
    version = ".".join(map(str, sys.version_info[:3]))
    if sys.version_info < (3, 10):  # noqa: UP036 - runtime guard for mis-installed interpreters
        return Check("python", FAIL, f"Python {version}", "AMR-Surveil needs Python 3.10 or newer.")
    if sys.version_info >= (3, 14):
        return Check(
            "python",
            WARN,
            f"Python {version}",
            "Binary wheels (pyarrow, numpy, pandas) can lag the newest Python. If an import fails, "
            "use a Python 3.11-3.13 environment: `conda create -n amr python=3.12`.",
        )
    return Check("python", OK, f"Python {version}")


def check_module(module: str, *, required: bool, purpose: str, fix: str) -> Check:
    version, error = _import_version(module)
    if version is not None:
        return Check(module, OK, version)
    return Check(module, FAIL if required else WARN, f"{purpose}: {error}", fix)


def check_smoke() -> Check:
    """Run the pipeline end to end on 24 synthetic isolates."""
    try:
        from amr_surveil.ui_support import run_demo

        with tempfile.TemporaryDirectory(prefix="amr_surveil_doctor_") as tmp:
            summary, out = run_demo(Path(tmp), n_samples=24)
            ok = (out / "report.md").is_file() and (out / "run_manifest.json").is_file()
            if not ok:
                return Check("pipeline smoke run", FAIL, "outputs missing after the run")
            return Check(
                "pipeline smoke run",
                OK,
                f"{summary.get('n_samples')} synthetic samples analysed",
            )
    except Exception as exc:
        return Check(
            "pipeline smoke run",
            FAIL,
            f"{type(exc).__name__}: {exc}",
            "Report this with the full traceback.",
        )


def run_checks(*, smoke: bool = True) -> list[Check]:
    checks = [
        check_python(),
        check_module("pandas", required=True, purpose="core tables", fix="pip install -U pandas"),
        check_module("numpy", required=True, purpose="numerics", fix="pip install -U numpy"),
        check_module("yaml", required=True, purpose="configuration", fix="pip install -U PyYAML"),
        check_module(
            "streamlit",
            required=False,
            purpose="web interface (the CLI works without it)",
            fix="pip install -U streamlit",
        ),
        check_module(
            "pyarrow",
            required=False,
            purpose="used by Streamlit to display tables; the web interface fails without it",
            fix="pip uninstall -y pyarrow && pip install -U pyarrow "
            "(on conda+pip mixes, prefer `conda install -c conda-forge pyarrow`)",
        ),
        check_module(
            "rpds",
            required=False,
            purpose="JSON-schema stack used by Streamlit (rpds-py); a broken wheel crashes some Streamlit elements",
            fix="pip uninstall -y rpds-py && pip install -U rpds-py "
            "(or use Python 3.11-3.13, where binary wheels are most mature)",
        ),
        check_module(
            "plotly",
            required=False,
            purpose="charts and maps in the web interface",
            fix="pip install -U plotly",
        ),
    ]
    if smoke:
        checks.append(check_smoke())
    return checks


def render(checks: list[Check]) -> str:
    lines = []
    for check in checks:
        lines.append(f"[{check.status:<4}] {check.name:<20} {check.detail}")
        if check.fix and check.status != OK:
            lines.append(f"         fix: {check.fix}")
    failed = sum(c.status == FAIL for c in checks)
    warned = sum(c.status == WARN for c in checks)
    lines.append("")
    lines.append(
        "All required checks passed." if not failed else f"{failed} required check(s) failed."
    )
    if warned:
        lines.append(f"{warned} optional component(s) need attention (see WARN lines).")
    return "\n".join(lines)


def exit_code(checks: list[Check]) -> int:
    return 1 if any(c.status == FAIL for c in checks) else 0
