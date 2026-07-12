#!/usr/bin/env python3
"""Step 2: select models on training data and evaluate held-out test files."""

from __future__ import annotations

import argparse
from pathlib import Path

import bootstrap  # noqa: F401
from stroke_rehab.evaluation import run_decoding


def run(data_dir: Path, output_dir: Path) -> None:
    results = run_decoding(data_dir, output_dir)
    print("\nSelected-model results")
    print(results["selected"].drop(columns=["confusion_matrix"]).to_string(index=False))
    print("\nFixed CSP+LDA results")
    print(results["fixed"].drop(columns=["confusion_matrix"]).to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    args = parser.parse_args()
    run(args.data_dir, args.output_dir)


if __name__ == "__main__":
    main()
