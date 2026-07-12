#!/usr/bin/env python3
"""Step 5: create PNG figures from the saved result tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import bootstrap  # noqa: F401
from stroke_rehab.plotting import make_all_python_figures


def run(output_dir: Path) -> None:
    make_all_python_figures(output_dir)
    print(f"Figures saved under: {output_dir / 'figures'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    args = parser.parse_args()
    run(args.output_dir)


if __name__ == "__main__":
    main()
