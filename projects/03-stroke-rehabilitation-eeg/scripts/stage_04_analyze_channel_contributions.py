#!/usr/bin/env python3
"""Step 4: estimate channel contributions with repeated trial shuffling."""

from __future__ import annotations

import argparse
from pathlib import Path

import bootstrap  # noqa: F401
from stroke_rehab.electrode_analysis import run_channel_contributions


def run(data_dir: Path, output_dir: Path, repeats: int = 30) -> None:
    results = run_channel_contributions(
        data_dir=data_dir,
        output_dir=output_dir,
        n_repeats=repeats,
    )
    print(results["contributions"].head(20).to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--repeats", type=int, default=30)
    args = parser.parse_args()
    run(args.data_dir, args.output_dir, repeats=args.repeats)


if __name__ == "__main__":
    main()
