"""Command-line interface for reproducible AMR-Surveil runs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml

from amr_surveil._version import __version__
from amr_surveil.adapters.amrfinder import normalize_amrfinder_reports
from amr_surveil.adapters.epiresnet import (
    normalize_epiresnet_predictions,
    normalize_epiresnet_test_predictions,
)
from amr_surveil.adapters.fishamr import (
    normalize_fishamr_calls,
    normalize_fishamr_feature_table,
    normalize_fishamr_plasmids,
    read_fishamr_warehouse,
)
from amr_surveil.doctor import exit_code, render, run_checks
from amr_surveil.exceptions import AMRSurveilError
from amr_surveil.pipeline import run_pipeline
from amr_surveil.provenance import verify_run, write_templates
from amr_surveil.synthetic import generate_synthetic_dataset

DEMO_CONFIG = {
    "analysis": {"deduplication_rule": "none", "reference_date": "2026-11-07"},
    "quality": {"strict": True},
    "study": {
        "name": "AMR-Surveil synthetic demonstration",
        "question": "Does the pipeline recover the planted gene-phenotype association and flag the planted method change?",
        "inclusion_criteria": "All synthetic isolates.",
        "exclusion_criteria": "None.",
        "sample_selection_documented": True,
        "selection_description": "Synthetic data generated with a fixed seed; no real isolates or sampling frame exist.",
        "predictions_independence_attested": False,
    },
}


def _add_inputs(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--samples", required=True, help="Sample/isolate manifest CSV")
    parser.add_argument("--amr-calls", help="Canonical AMR determinant-call CSV")
    parser.add_argument(
        "--screening",
        help="Screening provenance CSV (which samples were genotyped, with tool/database)",
    )
    parser.add_argument("--ast", help="Observed phenotypic AST CSV")
    parser.add_argument(
        "--breakpoints", help="Breakpoint table CSV used to re-check recorded S/I/R categories"
    )
    parser.add_argument("--aggregate", help="Separate aggregate surveillance counts CSV")
    parser.add_argument("--predictions", help="Separate model predictions CSV (e.g. EpiResNet)")
    parser.add_argument("--plasmids", help="Canonical plasmid/replicon results CSV")
    parser.add_argument(
        "--incidence", help="Separate population-based resistant-infection count CSV"
    )
    parser.add_argument("--consumption", help="Separate antimicrobial-consumption/DDD CSV")
    parser.add_argument("--gene-drug-map", help="Cited gene-to-antimicrobial mapping CSV")
    parser.add_argument("--config", help="YAML analysis/quality/study configuration")
    parser.add_argument("--out", required=True, help="Output directory")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="amr-surveil",
        description="Reproducible, provenance-aware AMR genomic and epidemiological surveillance workbench.",
    )
    parser.add_argument("--version", action="version", version=f"AMR-Surveil {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    _add_inputs(
        sub.add_parser(
            "run",
            help="Validate inputs and run the analysis (strict quality gate unless configured otherwise)",
        )
    )
    _add_inputs(
        sub.add_parser(
            "audit",
            help="Run in exploratory mode: emit every QA finding and outputs even if errors exist",
        )
    )
    _add_inputs(
        sub.add_parser(
            "qa", help="Exploratory run that exits with status 1 when ERROR findings exist (for CI)"
        )
    )

    demo = sub.add_parser(
        "demo", help="Generate deterministic synthetic inputs and run the full demo"
    )
    demo.add_argument("--out", required=True)
    demo.add_argument(
        "--samples", type=int, default=120, help="Number of synthetic samples (default: 120)"
    )
    demo.add_argument("--seed", type=int, default=20261107)
    demo.add_argument(
        "--inject-defects",
        action="store_true",
        help="Add documented data defects and run exploratory mode",
    )

    verify = sub.add_parser(
        "verify", help="Re-hash a finished run's outputs against its run_manifest.json"
    )
    verify.add_argument("results", help="Results directory containing run_manifest.json")
    verify.add_argument(
        "--check-inputs",
        action="store_true",
        help="Also re-hash the input files at the paths recorded in the manifest",
    )

    doctor = sub.add_parser(
        "doctor", help="Check the Python environment and run a tiny end-to-end smoke test"
    )
    doctor.add_argument(
        "--no-smoke", action="store_true", help="Only check imports; skip the smoke run"
    )

    public = sub.add_parser(
        "public-data",
        help="Retrieve or import openly available data as canonical AMR-Surveil inputs",
    )
    psub = public.add_subparsers(dest="public_command", required=True)
    live = psub.add_parser("bvbrc", help="Query BV-BRC live for isolate-level AST (bounded)")
    live.add_argument("--organism", required=True, help='Species, e.g. "Escherichia coli"')
    live.add_argument("--antibiotic", required=True, nargs="+", help="One or more antibiotic names")
    live.add_argument("--max-records", type=int, default=2000)
    live.add_argument("--year-from", type=int)
    live.add_argument("--year-to", type=int)
    live.add_argument("--country", nargs="+", help="Keep only these countries")
    live.add_argument("--cache-dir", default=None, help="Default: per-user cache folder")
    live.add_argument("--refresh", action="store_true", help="Ignore any cached copy")
    live.add_argument("--out", required=True)
    imp = psub.add_parser("import", help="Convert a file you downloaded from a public source")
    imp.add_argument("--source", required=True, choices=["bvbrc", "ncbi_pathogen", "ecdc_atlas"])
    imp.add_argument("--input", required=True, help="CSV or TSV downloaded from the source")
    imp.add_argument(
        "--include-genotypes",
        action="store_true",
        help="NCBI only: also import genotypes (exploratory)",
    )
    imp.add_argument("--out", required=True)
    psub.add_parser("sources", help="List supported sources with their terms")

    templates = sub.add_parser(
        "templates", help="Write header-only input templates and a default config.yaml"
    )
    templates.add_argument("--out", required=True)
    templates.add_argument("--force", action="store_true", help="Overwrite existing files")

    fish = sub.add_parser("normalize-fishamr", help="Normalize a FishAMR-Link export")
    fish.add_argument("--input", required=True)
    fish.add_argument("--output", required=True)
    fish.add_argument(
        "--kind", choices=["amr", "plasmid", "feature-table", "warehouse"], default="amr"
    )
    fish.add_argument("--mapping", help="JSON column mapping (kinds amr and plasmid)")
    fish.add_argument("--sample-id", help="Sample ID for a feature table that has none")
    fish.add_argument("--organism")
    fish.add_argument("--plasmids-output", help="Where to write replicon rows (feature-table)")
    fish.add_argument("--samples-output", help="Where to write the isolate manifest (warehouse)")
    fish.add_argument(
        "--include-seeded-baseline",
        action="store_true",
        help="Warehouse: keep the simulated baseline isolates",
    )

    amr = sub.add_parser(
        "normalize-amrfinder", help="Normalize AMRFinderPlus TSV report(s) to calls + screening"
    )
    amr.add_argument("--input", nargs="+", required=True)
    amr.add_argument("--output", required=True, help="AMR calls CSV")
    amr.add_argument("--screening-output", help="Screening CSV (needed for negative genotypes)")
    amr.add_argument(
        "--sample-id", nargs="*", help="Sample IDs, one per input, if the report has no Name column"
    )
    amr.add_argument("--organism")
    amr.add_argument(
        "--tool-version", help="AMRFinderPlus version, e.g. 3.12.8 (not stored in reports)"
    )
    amr.add_argument("--database-version", help="AMRFinderPlus database version, e.g. 2024-01-31.1")
    amr.add_argument(
        "--element-types", nargs="+", default=["AMR"], help="Element types to keep (default: AMR)"
    )

    epi = sub.add_parser("normalize-epiresnet", help="Normalize EpiResNet predictions")
    epi.add_argument("--input", required=True)
    epi.add_argument("--output", required=True)
    epi.add_argument("--mapping", help="JSON column mapping for a generic prediction CSV")
    epi.add_argument(
        "--test-predictions", action="store_true", help="Input is EpiResNet v5 test_predictions.csv"
    )
    epi.add_argument("--final-metrics", help="EpiResNet final_metrics.json (threshold/temperature)")
    epi.add_argument(
        "--run-manifest", help="EpiResNet run_manifest.json (split strategy, training-data hash)"
    )
    epi.add_argument("--checkpoint", help="Model checkpoint to fingerprint by SHA-256")
    return parser


def _run(args: argparse.Namespace, *, exploratory: bool = False) -> dict[str, object]:
    return run_pipeline(
        samples_path=args.samples,
        amr_calls_path=args.amr_calls,
        screening_path=args.screening,
        ast_path=args.ast,
        breakpoints_path=args.breakpoints,
        aggregate_path=args.aggregate,
        predictions_path=args.predictions,
        plasmids_path=args.plasmids,
        incidence_path=args.incidence,
        consumption_path=args.consumption,
        gene_drug_map_path=args.gene_drug_map,
        config_path=args.config,
        output_dir=args.out,
        config_overrides={"quality": {"strict": False}} if exploratory else None,
    )


def _short(summary: dict[str, object]) -> str:
    keys = [
        "n_samples",
        "n_ast_records",
        "n_amr_calls_used",
        "n_qa_errors",
        "n_qa_warnings",
        "genotype_phenotype_decision",
        "n_concordance_strata",
        "n_comparability_flags",
        "n_excluded_records",
    ]
    return json.dumps({k: summary.get(k) for k in keys}, indent=2)


def _public_data(args: argparse.Namespace) -> int:
    from datetime import datetime, timezone

    from amr_surveil.datasources import (
        SOURCES,
        default_cache_dir,
        fetch_bvbrc,
        parse_bvbrc_export,
        parse_ecdc_atlas,
        parse_pathogen_table,
        write_inputs,
    )
    from amr_surveil.datasources.common import read_table_upload

    if args.public_command == "sources":
        for info in SOURCES.values():
            print(
                f"{info.source_id} [{info.access}] {info.name}\n  terms ({info.terms_checked}): {info.terms}\n  {info.terms_url}"
            )
        return 0
    if args.public_command == "bvbrc":
        dataset = fetch_bvbrc(
            args.organism,
            args.antibiotic,
            cache_dir=args.cache_dir or default_cache_dir(),
            refresh=args.refresh,
            max_records=args.max_records,
            year_from=args.year_from,
            year_to=args.year_to,
            countries=args.country,
        )
    else:
        path = Path(args.input)
        table = read_table_upload(path.read_bytes(), path.name)
        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d")
        if args.source == "bvbrc":
            dataset = parse_bvbrc_export(table, retrieved_at=stamp)
        elif args.source == "ncbi_pathogen":
            dataset = parse_pathogen_table(
                table, include_genotypes=args.include_genotypes, retrieved_at=stamp
            )
        else:
            dataset = parse_ecdc_atlas(table, retrieved_at=stamp)
    paths = write_inputs(dataset, args.out)
    print(json.dumps(dataset.summary(), indent=2))
    for key, written in paths.items():
        print(f"{key}: {written}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "run":
            summary = _run(args)
            print(_short(summary))
            print(f"\nCompleted. Report: {Path(args.out) / 'report.md'}")
        elif args.command in {"audit", "qa"}:
            summary = _run(args, exploratory=True)
            print(_short(summary))
            print(
                f"\nExploratory run completed. QA findings: {Path(args.out) / 'audit_issues.csv'}"
            )
            if args.command == "qa" and summary["n_qa_errors"]:
                return 1
        elif args.command == "demo":
            root = Path(args.out)
            inputs = generate_synthetic_dataset(
                root / "inputs",
                n_samples=args.samples,
                seed=args.seed,
                inject_defects=args.inject_defects,
            )
            config = root / "inputs" / "config.yaml"
            config.write_text(yaml.safe_dump(DEMO_CONFIG, sort_keys=False), encoding="utf-8")
            summary = run_pipeline(
                samples_path=inputs["samples"],
                amr_calls_path=inputs["amr_calls"],
                screening_path=inputs["screening"],
                ast_path=inputs["ast"],
                breakpoints_path=inputs["breakpoints"],
                aggregate_path=inputs["aggregate"],
                predictions_path=inputs["predictions"],
                plasmids_path=inputs["plasmids"],
                incidence_path=inputs["incidence"],
                consumption_path=inputs["consumption"],
                gene_drug_map_path=inputs["gene_drug_map"],
                config_path=config,
                output_dir=root / "results",
                config_overrides={"quality": {"strict": False}} if args.inject_defects else None,
            )
            print(_short(summary))
            print(f"\nSynthetic demo report: {root / 'results' / 'report.md'}")
        elif args.command == "verify":
            report = verify_run(args.results, check_inputs=args.check_inputs)
            print(json.dumps(report, indent=2))
            if not report["ok"]:
                return 1
        elif args.command == "public-data":
            return _public_data(args)
        elif args.command == "doctor":
            checks = run_checks(smoke=not args.no_smoke)
            print(render(checks))
            return exit_code(checks)
        elif args.command == "templates":
            for path in write_templates(args.out, force=args.force):
                print(path)
        elif args.command == "normalize-fishamr":
            if args.kind in {"amr", "plasmid"}:
                if not args.mapping:
                    raise AMRSurveilError("--mapping is required for --kind amr/plasmid")
                fn = normalize_fishamr_calls if args.kind == "amr" else normalize_fishamr_plasmids
                print(
                    f"Normalized {len(fn(args.input, args.output, args.mapping))} FishAMR-Link row(s) to {args.output}"
                )
            elif args.kind == "feature-table":
                print(
                    json.dumps(
                        normalize_fishamr_feature_table(
                            args.input,
                            args.output,
                            args.plasmids_output,
                            sample_id=args.sample_id,
                            organism=args.organism,
                        ),
                        indent=2,
                    )
                )
            else:
                if not args.samples_output:
                    raise AMRSurveilError("--samples-output is required for --kind warehouse")
                print(
                    json.dumps(
                        read_fishamr_warehouse(
                            args.input,
                            args.samples_output,
                            args.output,
                            include_seeded_baseline=args.include_seeded_baseline,
                        ),
                        indent=2,
                    )
                )
        elif args.command == "normalize-amrfinder":
            print(
                json.dumps(
                    normalize_amrfinder_reports(
                        args.input,
                        args.output,
                        args.screening_output,
                        sample_ids=args.sample_id or None,
                        organism=args.organism,
                        tool_version=args.tool_version,
                        database_version=args.database_version,
                        element_types=tuple(args.element_types),
                    ),
                    indent=2,
                )
            )
        elif args.command == "normalize-epiresnet":
            if args.test_predictions:
                print(
                    json.dumps(
                        normalize_epiresnet_test_predictions(
                            args.input,
                            args.output,
                            final_metrics_path=args.final_metrics,
                            run_manifest_path=args.run_manifest,
                            checkpoint_path=args.checkpoint,
                        ),
                        indent=2,
                    )
                )
            else:
                if not args.mapping:
                    raise AMRSurveilError(
                        "--mapping is required unless --test-predictions is given"
                    )
                frame = normalize_epiresnet_predictions(args.input, args.output, args.mapping)
                print(f"Normalized {len(frame)} EpiResNet prediction row(s) to {args.output}")
        return 0
    except (AMRSurveilError, ValueError, OSError, yaml.YAMLError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
