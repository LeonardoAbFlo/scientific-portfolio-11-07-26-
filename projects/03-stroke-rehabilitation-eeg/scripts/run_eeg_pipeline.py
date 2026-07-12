#!/usr/bin/env python3
"""Run the complete EEG decoding pipeline in sequence."""

from __future__ import annotations

import argparse
from pathlib import Path

from stage_01_validate_eeg_data import run as validate
from stage_02_run_decoding_analysis import run as decode
from stage_03_export_diagnostics import run as diagnostics
from stage_04_analyze_channel_contributions import run as channel_contributions
from stage_05_create_figures import run as figures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument(
        "--skip-channel-contributions",
        action="store_true",
        help="Skip the slower permutation-based electrode analysis.",
    )
    args = parser.parse_args()

    validate(args.data_dir, args.output_dir, strict=True)
    decode(args.data_dir, args.output_dir)
    diagnostics(args.data_dir, args.output_dir)
    if not args.skip_channel_contributions:
        channel_contributions(args.data_dir, args.output_dir, repeats=args.repeats)
    figures(args.output_dir)
    print(f"\nPipeline complete. Outputs: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
