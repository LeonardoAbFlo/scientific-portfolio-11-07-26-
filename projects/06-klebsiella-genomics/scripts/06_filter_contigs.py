#!/usr/bin/env python3
"""Create a filtered contig set and an auditable keep/exclude decision table."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from collections import defaultdict
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Iterable


SAMPLE_COLUMNS = ("sample_id", "barcode", "UCL_ID", "Sample")
CONTIG_COLUMNS = ("contig_id", "Contig", "Name")
MOLECULE_COLUMNS = ("molecule_type", "Molecule_type")
CLASSIFICATION_COLUMNS = ("classification", "Classification")
NEIGHBOR_COLUMNS = ("mash_neighbor_identification", "Mash_neighbor")
CIRCULAR_COLUMNS = ("circ", "circularity_status", "Circular")


@dataclass(frozen=True)
class Decision:
    sample_id: str
    contig_id: str
    molecule_type: str
    expected_genus: str
    detected_genus: str
    mixed_sample: bool
    decision: str
    reason: str


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        return reader.fieldnames, [dict(row) for row in reader]


def select_column(headers: Iterable[str], candidates: Iterable[str], label: str) -> str:
    header_set = set(headers)
    for candidate in candidates:
        if candidate in header_set:
            return candidate
    raise ValueError(f"Cannot find {label} column; tried {', '.join(candidates)}")


def genus(value: str | None) -> str:
    if not value:
        return ""
    raw_text = value.strip()
    gtdb = re.search(r"(?:^|;)g__([A-Za-z][A-Za-z-]*)", raw_text)
    if gtdb:
        return gtdb.group(1)
    text = raw_text.replace("_", " ")
    if text.lower() in {"unknown", "unclassified", "na", "n/a"}:
        return ""
    first = re.search(r"[A-Za-z][A-Za-z-]*", text)
    return first.group(0) if first else ""


def safe_identifier(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")


def similar(left: str, right: str, threshold: float) -> bool:
    if not left or not right:
        return False
    return SequenceMatcher(None, left.casefold(), right.casefold()).ratio() >= threshold


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--taxonomy", required=True, type=Path)
    parser.add_argument("--contigs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--decisions", type=Path)
    parser.add_argument("--fuzzy-threshold", type=float, default=0.70)
    parser.add_argument("--exclude-cross-genus-plasmids", action="store_true")
    parser.add_argument("--exclude-noncircular-kp", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not 0 <= args.fuzzy_threshold <= 1:
        raise SystemExit("--fuzzy-threshold must be between 0 and 1")
    metadata_headers, metadata_rows = read_tsv(args.metadata)
    taxonomy_headers, taxonomy_rows = read_tsv(args.taxonomy)

    metadata_sample = select_column(metadata_headers, SAMPLE_COLUMNS, "metadata sample")
    species_column = select_column(
        metadata_headers,
        ("expected_species", "Species_Greece", "species", "Species"),
        "expected species",
    )
    sample_column = select_column(taxonomy_headers, SAMPLE_COLUMNS, "taxonomy sample")
    contig_column = select_column(taxonomy_headers, CONTIG_COLUMNS, "contig")
    molecule_column = select_column(taxonomy_headers, MOLECULE_COLUMNS, "molecule type")
    classification_column = select_column(taxonomy_headers, CLASSIFICATION_COLUMNS, "classification")
    neighbor_column = next((c for c in NEIGHBOR_COLUMNS if c in taxonomy_headers), "")
    circular_column = next((c for c in CIRCULAR_COLUMNS if c in taxonomy_headers), "")

    expected: dict[str, str] = {}
    for row in metadata_rows:
        sample_id = row[metadata_sample].strip()
        if sample_id in expected:
            raise ValueError(f"Duplicate metadata sample: {sample_id}")
        expected[sample_id] = genus(row[species_column])

    chromosome_genera: dict[str, set[str]] = defaultdict(set)
    for row in taxonomy_rows:
        if row[molecule_column].strip().casefold() == "chromosome":
            detected = genus(row[classification_column])
            if detected:
                chromosome_genera[row[sample_column].strip()].add(detected)

    decisions: list[Decision] = []
    kept = excluded = 0
    args.output.mkdir(parents=True, exist_ok=True)

    for row in taxonomy_rows:
        sample_id = row[sample_column].strip()
        contig_id = row[contig_column].strip()
        molecule = row[molecule_column].strip().casefold()
        observed_text = row[classification_column] if molecule == "chromosome" else row.get(neighbor_column, "")
        detected = genus(observed_text)
        expected_genus = expected.get(sample_id, "")
        mixed = len(chromosome_genera.get(sample_id, set())) > 1
        decision = "keep"
        reason = "classification consistent or insufficient evidence to exclude"

        if molecule == "chromosome" and expected_genus and detected:
            if not similar(detected, expected_genus, args.fuzzy_threshold):
                decision = "exclude"
                reason = "chromosome genus conflicts with expected sample genus"
        if (
            decision == "keep"
            and args.exclude_noncircular_kp
            and molecule == "chromosome"
            and re.search(r"Klebsiella[ _]pneumoniae", row[classification_column], re.I)
            and circular_column
            and row.get(circular_column, "").strip().upper() not in {"Y", "YES", "CIRCULAR"}
        ):
            decision = "exclude"
            reason = "non-circular K. pneumoniae chromosome excluded by requested rule"
        if (
            decision == "keep"
            and args.exclude_cross_genus_plasmids
            and molecule == "plasmid"
            and mixed
            and expected_genus
            and detected
            and not similar(detected, expected_genus, args.fuzzy_threshold)
        ):
            decision = "exclude"
            reason = "cross-genus plasmid in a mixed sample excluded by requested rule"

        decisions.append(
            Decision(
                sample_id,
                contig_id,
                molecule,
                expected_genus,
                detected,
                mixed,
                decision,
                reason,
            )
        )

        if decision == "exclude":
            excluded += 1
            continue
        source = args.contigs / safe_identifier(sample_id) / f"{safe_identifier(contig_id)}.fasta"
        if not source.is_file():
            raise FileNotFoundError(f"Retained contig FASTA not found: {source}")
        destination = args.output / safe_identifier(sample_id) / source.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and not args.force:
            raise FileExistsError(f"Refusing to overwrite {destination}; use --force")
        shutil.copy2(source, destination)
        kept += 1

    decision_path = args.decisions or args.output / "filter_decisions.tsv"
    decision_path.parent.mkdir(parents=True, exist_ok=True)
    with decision_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(Decision.__annotations__), delimiter="\t")
        writer.writeheader()
        writer.writerows(decision.__dict__ for decision in decisions)

    print(f"Kept {kept} contigs; excluded {excluded}; wrote {decision_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
