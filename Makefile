.PHONY: install test lint format check demo demo-defects schemas benchmark docker

install:
	python -m pip install -e ".[dev]"

test:
	pytest --cov=amr_surveil --cov-report=term-missing

lint:
	ruff check .
	ruff format --check .

format:
	ruff format .

check: lint
	python -m compileall -q amr_surveil app.py benchmarks scripts
	pytest -q

demo:
	amr-surveil demo --out results/demo --samples 120 --seed 20261107

demo-defects:
	amr-surveil demo --out results/demo_defects --samples 40 --inject-defects

schemas:
	python scripts/generate_schemas.py

benchmark:
	python benchmarks/run_benchmarks.py --sizes 100 1000 5000 20000 --output benchmarks/results/latest.json

docker:
	docker build -t amr-surveil:1.1.1 .
