# Contributing

Contributions that improve reproducibility, scientific transparency, data validation, accessibility or test coverage are welcome.

## Before opening a pull request

1. Open an issue for significant changes to data contracts or analysis behaviour.
2. Keep the default workflow conservative: no fuzzy isolate joins, no inferred negative gene calls, no clinical interpretations.
3. Add or update unit tests for every behaviour change.
4. Use synthetic fixtures only; do not submit patient-identifiable, confidential or restricted data.
5. Record changes to methods, outputs, defaults and compatibility in `CHANGELOG.md`.

## Local development

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
python -m pip install -e ".[dev]"
pytest
ruff check .
python -m compileall -q amr_surveil app.py
```

## Scientific review expectations

Any new organism-specific method, breakpoint interpretation, gene–drug mapping, deduplication rule, genomic distance or One Health comparison must include authoritative references, scope, versioning, limitations and test fixtures. A method should not be described as clinically validated without a suitable independent validation design and evidence.
