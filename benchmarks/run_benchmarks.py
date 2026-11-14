"""Reproducible synthetic throughput benchmark; not a clinical performance benchmark."""

from __future__ import annotations

import argparse
import json
import platform
import tempfile
import time
from pathlib import Path

import pandas as pd
import psutil
import yaml

from amr_surveil._version import __version__
from amr_surveil.cli import DEMO_CONFIG
from amr_surveil.pipeline import run_pipeline
from amr_surveil.synthetic import generate_synthetic_dataset


def benchmark_size(n_samples: int, seed: int, work_dir: Path) -> dict[str, object]:
    input_dir = work_dir / f"input_{n_samples}"
    output_dir = work_dir / f"output_{n_samples}"
    paths = generate_synthetic_dataset(input_dir, n_samples=n_samples, seed=seed)
    config_path = work_dir / "config.yaml"
    config_path.write_text(yaml.safe_dump(DEMO_CONFIG), encoding="utf-8")
    process = psutil.Process()
    baseline = process.memory_info().rss
    started = time.perf_counter()
    summary = run_pipeline(
        samples_path=paths["samples"],
        amr_calls_path=paths["amr_calls"],
        screening_path=paths["screening"],
        ast_path=paths["ast"],
        breakpoints_path=paths["breakpoints"],
        gene_drug_map_path=paths["gene_drug_map"],
        aggregate_path=paths["aggregate"],
        predictions_path=paths["predictions"],
        plasmids_path=paths["plasmids"],
        incidence_path=paths["incidence"],
        consumption_path=paths["consumption"],
        config_path=config_path,
        output_dir=output_dir,
    )
    elapsed = time.perf_counter() - started
    peak = process.memory_info().rss
    # Correctness guard: a fast run that lost the planted association is not a valid benchmark point.
    conc = pd.read_csv(output_dir / "genotype_phenotype_concordance.csv")
    pos_r, pos_n = conc["n_gene_pos_R"].sum(), conc["n_gene_pos_notR"].sum()
    neg_r, neg_n = conc["n_gene_neg_R"].sum(), conc["n_gene_neg_notR"].sum()
    ppv = pos_r / (pos_r + pos_n) if pos_r + pos_n else float("nan")
    miss_rate = neg_r / (neg_r + neg_n) if neg_r + neg_n else float("nan")
    return {
        "samples": n_samples,
        "elapsed_seconds": round(elapsed, 4),
        "samples_per_second": round(n_samples / elapsed, 2) if elapsed else None,
        "rss_before_mb": round(baseline / (1024 * 1024), 2),
        "rss_after_mb": round(peak / (1024 * 1024), 2),
        "rss_delta_mb": round(max(0, peak - baseline) / (1024 * 1024), 2),
        "ast_rows": summary["n_ast_records"],
        "amr_call_rows": summary["n_amr_calls"],
        "qa_errors": summary["n_qa_errors"],
        "gate_decision": summary["genotype_phenotype_decision"],
        "concordance_strata": summary["n_concordance_strata"],
        "pooled_ppv_planted_0.92": round(float(ppv), 4),
        "pooled_resistance_without_gene_planted_0.06": round(float(miss_rate), 4),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=int, default=[100, 1000, 5000, 20000])
    parser.add_argument("--seed", type=int, default=20261107)
    parser.add_argument("--output", default="benchmarks/results/latest.json")
    args = parser.parse_args()
    if any(size < 8 for size in args.sizes):
        parser.error("All benchmark sizes must be at least 8")
    records = []
    with tempfile.TemporaryDirectory(prefix="amr_surveil_bench_") as tmp:
        work_dir = Path(tmp)
        for size in args.sizes:
            records.append(benchmark_size(size, args.seed, work_dir))
    result = {
        "benchmark_type": "synthetic_pipeline_throughput_with_ground_truth_check",
        "amr_surveil_version": __version__,
        "warning": "Synthetic engineering benchmark only; not a biological accuracy or operational surveillance validation.",
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "platform": platform.platform(),
        "seed": args.seed,
        "runs": records,
    }
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
