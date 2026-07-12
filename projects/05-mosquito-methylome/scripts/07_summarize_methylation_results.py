#!/usr/bin/env python3
"""Build coverage-aware sample and top-site summaries from Modkit BEDMethyl."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Site:
    chrom: str
    position_1based: int
    strand: str
    coverage: int
    modified_reads: int
    fraction_modified: float


def parse_int(value: str) -> int:
    return int(float(value))


def read_sites(path: Path) -> tuple[list[Site], int]:
    sites: list[Site] = []
    malformed = 0
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            try:
                if len(row) < 12:
                    raise ValueError("BEDMethyl row has fewer than 12 columns")
                start = parse_int(row[1])
                coverage = parse_int(row[9])
                percent_modified = float(row[10])
                modified_reads = parse_int(row[11])
                sites.append(
                    Site(
                        chrom=row[0],
                        position_1based=start + 1,
                        strand=row[5],
                        coverage=coverage,
                        modified_reads=modified_reads,
                        fraction_modified=percent_modified / 100.0,
                    )
                )
            except (IndexError, TypeError, ValueError):
                malformed += 1
    return sites, malformed


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "species", "barcode", "reference_fasta"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Sample sheet must contain: {', '.join(sorted(required))}")
        rows = [dict(row) for row in reader if row.get("sample_id") and not row["sample_id"].startswith("#")]
    sample_ids = [row["sample_id"] for row in rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Sample IDs must be unique")
    return rows


def motif_summary(path: Path) -> tuple[str, int, str]:
    if not path.exists():
        return "missing", 0, ""
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if row]
    if len(rows) <= 1:
        return "not_detected", 0, ""
    header = rows[0]
    data = rows[1:]
    motif_index = next((header.index(name) for name in ("motif", "Motif", "sequence") if name in header), 0)
    return "detected", len(data), data[0][motif_index]


def fmt(value: float | int | None) -> str | int:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    if isinstance(value, float):
        return f"{value:.6f}"
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--pileups", required=True, type=Path)
    parser.add_argument("--motifs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-coverage", type=int, default=20)
    parser.add_argument("--candidate-fraction", type=float, default=0.20)
    parser.add_argument("--weak-fraction", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=10)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.min_coverage < 1 or args.top_n < 1:
        raise SystemExit("Coverage and top-n values must be positive")
    if not 0 <= args.weak_fraction <= args.candidate_fraction <= 1:
        raise SystemExit("Require 0 <= weak fraction <= candidate fraction <= 1")

    args.output.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, object]] = []
    top_rows: list[dict[str, object]] = []

    for sample in read_samples(args.samples):
        sample_id = sample["sample_id"]
        pileup = args.pileups / f"{sample_id}.bedmethyl"
        motif = args.motifs / f"{sample_id}.tsv"
        if pileup.exists():
            sites, malformed = read_sites(pileup)
            status = "ok" if sites else "empty"
        else:
            sites, malformed, status = [], 0, "missing"

        adequate = [site for site in sites if site.coverage >= args.min_coverage]
        candidates = [site for site in adequate if site.fraction_modified >= args.candidate_fraction]
        weak = [site for site in adequate if site.fraction_modified >= args.weak_fraction]
        ranked = sorted(adequate, key=lambda site: (site.fraction_modified, site.modified_reads, site.coverage), reverse=True)
        top = ranked[0] if ranked else None
        motif_status, motif_rows, top_motif = motif_summary(motif)

        if status == "missing":
            interpretation = "pileup_missing"
        elif not adequate:
            interpretation = "insufficient_coverage"
        elif candidates:
            interpretation = "candidate_modified_sites_detected"
        elif weak:
            interpretation = "weak_modified_site_signal"
        else:
            interpretation = "no_site_above_threshold"

        summary_rows.append(
            {
                "sample_id": sample_id,
                "species": sample["species"],
                "barcode": sample["barcode"],
                "pileup_status": status,
                "positions_evaluated": len(sites),
                "positions_adequate_coverage": len(adequate),
                "total_valid_coverage": sum(site.coverage for site in sites),
                "total_modified_reads": sum(site.modified_reads for site in sites),
                "mean_fraction_modified": fmt(sum(site.fraction_modified for site in sites) / len(sites) if sites else None),
                "max_fraction_modified_adequate_coverage": fmt(top.fraction_modified if top else None),
                "candidate_sites": len(candidates),
                "weak_sites": len(weak),
                "top_site": f"{top.chrom}:{top.position_1based}:{top.strand}" if top else "NA",
                "top_site_coverage": fmt(top.coverage if top else None),
                "top_site_modified_reads": fmt(top.modified_reads if top else None),
                "motif_status": motif_status,
                "motif_rows": motif_rows,
                "top_motif": top_motif or "NA",
                "malformed_pileup_rows": malformed,
                "interpretation": interpretation,
                "pileup_source": str(pileup),
                "motif_source": str(motif),
            }
        )

        for rank, site in enumerate(ranked[: args.top_n], start=1):
            top_rows.append(
                {
                    "sample_id": sample_id,
                    "species": sample["species"],
                    "rank": rank,
                    "chrom": site.chrom,
                    "position_1based": site.position_1based,
                    "strand": site.strand,
                    "coverage": site.coverage,
                    "modified_reads": site.modified_reads,
                    "fraction_modified": fmt(site.fraction_modified),
                }
            )

    for filename, rows in (("methylation_summary.tsv", summary_rows), ("top_sites.tsv", top_rows)):
        path = args.output / filename
        fieldnames = list(rows[0]) if rows else (["sample_id"] if filename == "top_sites.tsv" else [])
        with path.open("w", encoding="utf-8", newline="") as handle:
            if fieldnames:
                writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
    print(f"Wrote summaries for {len(summary_rows)} samples to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

