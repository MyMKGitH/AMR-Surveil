"""Documentation must not drift from the code."""

import re
from pathlib import Path

from amr_surveil import __version__
from amr_surveil import constants as C

ROOT = Path(__file__).resolve().parents[1]


def test_every_contract_column_is_documented():
    text = (ROOT / "docs" / "DATA_DICTIONARY.md").read_text(encoding="utf-8")
    names = (
        "SAMPLE_COLUMNS AMR_CALL_COLUMNS SCREENING_COLUMNS AST_COLUMNS BREAKPOINT_COLUMNS AGGREGATE_COLUMNS "
        "INCIDENCE_COLUMNS CONSUMPTION_COLUMNS PREDICTION_COLUMNS PLASMID_COLUMNS GENE_MAP_COLUMNS"
    ).split()
    missing = [c for n in names for c in getattr(C, n) if f"`{c}`" not in text]
    assert not missing, f"undocumented columns: {sorted(set(missing))}"


def test_every_pipeline_output_is_documented():
    from amr_surveil.pipeline import OUTPUT_FILES

    text = (ROOT / "docs" / "DATA_DICTIONARY.md").read_text(encoding="utf-8")
    assert not [f for f in OUTPUT_FILES if f != "report.md" and f"`{f}`" not in text]


def test_version_is_consistent_everywhere():
    assert __version__ == "1.1.1"
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert re.search(rf'^version = "{re.escape(__version__)}"', pyproject, re.M)
    assert f"version: {__version__}" in (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert f"## [{__version__}]" in (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert f"amr-surveil:{__version__}" in (ROOT / "docker-compose.yml").read_text(encoding="utf-8")
    assert f"ARG VERSION={__version__}" in (ROOT / "Dockerfile").read_text(encoding="utf-8")


def _documented_commands():
    import shlex

    commands = []
    for path in [ROOT / "README.md", *sorted((ROOT / "docs").glob("*.md"))]:
        for block in re.findall(
            r"```(?:bash|sh|console)?\n(.*?)```", path.read_text("utf-8"), re.S
        ):
            joined = re.sub(r"\\\n\s*", " ", block)
            for line in joined.splitlines():
                line = line.split("#")[0].strip()
                if line.startswith("amr-surveil ") and "--help" not in line:
                    commands.append((path.name, shlex.split(line)[1:]))
    return commands


def test_every_documented_cli_command_parses():
    from amr_surveil.cli import build_parser

    commands = _documented_commands()
    assert len(commands) >= 6
    parser = build_parser()
    for source, argv in commands:
        try:
            parser.parse_args(argv)
        except SystemExit as exc:  # argparse exits on error
            raise AssertionError(f"{source}: documented command does not parse: {argv}") from exc
