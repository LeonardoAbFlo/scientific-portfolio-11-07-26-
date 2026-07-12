#!/usr/bin/env python3
"""Step 3: export trigger and electrode-level diagnostic tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import bootstrap  # noqa: F401
from stroke_rehab.diagnostics import run_diagnostics


def run(data_dir: Path, output_dir: Path) -> None:
    results = run_diagnostics(data_dir, output_dir)
    for name, frame in results.items():
        print(f"{name}: {len(frame):,} rows")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    args = parser.parse_args()
    run(args.data_dir, args.output_dir)


if __name__ == "__main__":
    main()
