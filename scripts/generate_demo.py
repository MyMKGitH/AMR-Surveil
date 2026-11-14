"""Generate the bundled synthetic example dataset (and optionally a defective variant)."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from amr_surveil.cli import DEMO_CONFIG
from amr_surveil.synthetic import generate_synthetic_dataset


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=root / "data" / "demo")
    parser.add_argument("--samples", type=int, default=120)
    parser.add_argument("--seed", type=int, default=20261107)
    parser.add_argument(
        "--inject-defects", action="store_true", help="add the documented data defects"
    )
    args = parser.parse_args()
    paths = generate_synthetic_dataset(
        args.out, n_samples=args.samples, seed=args.seed, inject_defects=args.inject_defects
    )
    if not args.inject_defects:
        (args.out / "config.yaml").write_text(
            yaml.safe_dump(DEMO_CONFIG, sort_keys=False), encoding="utf-8"
        )
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
