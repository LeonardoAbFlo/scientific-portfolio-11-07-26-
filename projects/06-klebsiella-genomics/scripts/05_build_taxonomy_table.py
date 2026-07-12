#!/usr/bin/env python3
"""Combine contig manifest, MOB-suite, and GTDB-Tk results into one TSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def first_value(row: dict[str, str], *names: str) -> str:
    for name in names:
        if name in row and row[name] is not None:
            return row[name]
    return ""


def clean_contig(value: str) -> str:
    return Path(value.strip()).stem


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--mobsuite", required=True, type=Path)
    parser.add_argument("--gtdbtk", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = read_tsv(args.manifest)

    mob: dict[tuple[str, str], dict[str, str]] = {}
    for report in sorted(args.mobsuite.glob("*/contig_report.txt")):
        sample = report.parent.name
        for row in read_tsv(report):
            contig = clean_contig(first_value(row, "contig_id", "Contig", "contig"))
            if contig:
                mob[(sample, contig)] = row

    gtdb: dict[tuple[str, str], str] = {}
    for report in sorted(args.gtdbtk.glob("*/**/gtdbtk.*.summary.tsv")):
        sample = report.parts[len(args.gtdbtk.parts)]
        for row in read_tsv(report):
            contig = clean_contig(first_value(row, "user_genome", "Genome", "Name"))
            classification = first_value(row, "classification", "Classification")
            if contig:
                gtdb[(sample, contig)] = classification

    columns = [
        "sample_id",
        "contig_id",
        "length_bp",
        "molecule_type",
        "circularity_status",
        "mash_neighbor_identification",
        "rep_type(s)",
        "classification",
        "contig_fasta",
    ]
    output_rows: list[dict[str, str]] = []
    for row in manifest:
        sample = row["sample_id"]
        contig = clean_contig(row["contig_id"])
        mob_row = mob.get((sample, contig), {})
        output_rows.append(
            {
                "sample_id": sample,
                "contig_id": row["contig_id"],
                "length_bp": row["length_bp"],
                "molecule_type": first_value(mob_row, "molecule_type", "Molecule_type"),
                "circularity_status": first_value(mob_row, "circularity_status", "Circular"),
                "mash_neighbor_identification": first_value(mob_row, "mash_neighbor_identification"),
                "rep_type(s)": first_value(mob_row, "rep_type(s)", "rep_type"),
                "classification": gtdb.get((sample, contig), ""),
                "contig_fasta": row["contig_fasta"],
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {len(output_rows)} contigs to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

