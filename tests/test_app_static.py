"""Static guards for the Streamlit script.

A real Streamlit runtime is not available in CI for this suite, so these tests encode the two defects
a user hit on first local run, and would have caught both:

* Streamlit "magic" wraps any bare expression statement in ``st.write``. A conditional expression
  used as a statement (``st.dataframe(x) if cond else st.success(y)``) is therefore re-parsed from
  its source text and raised ``SyntaxError: '(' was never closed``.
* ``use_container_width`` is deprecated in current Streamlit.
"""

from __future__ import annotations

import ast
from pathlib import Path

APP = Path(__file__).resolve().parents[1] / "app.py"
UI_FILES = [APP, *sorted((APP.parent / "amr_surveil").glob("ui_*.py"))]


def _tree(path: Path) -> ast.AST:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def test_no_bare_expression_statements_that_trigger_streamlit_magic() -> None:
    offenders = []
    for path in UI_FILES:
        for node in ast.walk(_tree(path)):
            if isinstance(node, ast.Expr) and isinstance(
                node.value, (ast.IfExp, ast.BoolOp, ast.Compare, ast.BinOp, ast.Name, ast.Attribute)
            ):
                offenders.append(f"{path.name}:{node.lineno}")
    assert not offenders, f"bare expression statements (Streamlit magic) at {offenders}"


def test_deprecated_container_width_only_used_in_the_compat_fallback() -> None:
    offenders = []
    for path in UI_FILES:
        tree = _tree(path)
        allowed = {
            n.lineno
            for n in ast.walk(tree)
            if isinstance(n, ast.FunctionDef) and n.name == "_wide"
            for n in ast.walk(n)
            if hasattr(n, "lineno")
        }
        for node in ast.walk(tree):
            if isinstance(node, ast.keyword) and node.arg == "use_container_width":
                if node.value.lineno not in allowed:
                    offenders.append(f"{path.name}:{node.value.lineno}")
    assert not offenders, f"use_container_width outside the _wide fallback at {offenders}"


def test_app_script_compiles() -> None:
    compile(APP.read_text(encoding="utf-8"), str(APP), "exec")


def test_app_avoids_altair_backed_streamlit_charts():
    """st.bar_chart/line_chart/area_chart import altair -> jsonschema -> rpds-py; a broken
    rpds-py wheel (seen on Python 3.14) then crashes the tab. Use Plotly instead."""
    from pathlib import Path

    source = (Path(__file__).resolve().parents[1] / "app.py").read_text(encoding="utf-8")
    for name in (
        "st.bar_chart",
        "st.line_chart",
        "st.area_chart",
        "st.scatter_chart",
        "st.altair_chart",
    ):
        assert name not in source
