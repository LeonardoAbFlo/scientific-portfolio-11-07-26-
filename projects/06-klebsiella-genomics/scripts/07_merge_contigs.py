#!/usr/bin/env python3
"""Merge retained per-contig FASTA files into one assembly per sample."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.contigs.is_dir():
        raise SystemExit(f"Contig directory not found: {args.contigs}")
    args.output.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []

    for sample_dir in sorted(path for path in args.contigs.iterdir() if path.is_dir()):
        fasta_files = sorted(sample_dir.glob("*.fasta"))
        if not fasta_files:
            print(f"[WARN] No retained contigs for {sample_dir.name}")
            continue
        destination = args.output / sample_dir.name / "assembly.fasta"
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and not args.force:
            raise FileExistsError(f"Refusing to overwrite {destination}; use --force")
        total_bp = 0
        with destination.open("w", encoding="utf-8") as output:
            for fasta in fasta_files:
                text = fasta.read_text(encoding="utf-8").strip()
                if not text.startswith(">"):
                    raise ValueError(f"Invalid FASTA file: {fasta}")
                output.write(text + "\n")
                total_bp += sum(len(line.strip()) for line in text.splitlines() if not line.startswith(">"))
        manifest_rows.append(
            {
                "sample_id": sample_dir.name,
                "retained_contigs": len(fasta_files),
                "assembly_length_bp": total_bp,
                "assembly_fasta": str(destination.resolve()),
            }
        )

    manifest = args.output / "filtered_assembly_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "retained_contigs", "assembly_length_bp", "assembly_fasta"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Wrote {len(manifest_rows)} filtered assemblies and {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
