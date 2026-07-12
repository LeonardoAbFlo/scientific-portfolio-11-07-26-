#!/usr/bin/env python3
"""Step 1: validate MATLAB files, EEG dimensions and trigger labels."""

from __future__ import annotations

import argparse
from pathlib import Path

import bootstrap  # noqa: F401
from stroke_rehab.validation import validate_dataset


def run(data_dir: Path, output_dir: Path, strict: bool = True) -> None:
    report = validate_dataset(data_dir, output_dir, strict=strict)
    print(report.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Write the validation report without failing on missing/invalid files.",
    )
    args = parser.parse_args()
    run(args.data_dir, args.output_dir, strict=not args.allow_missing)


if __name__ == "__main__":
    main()
