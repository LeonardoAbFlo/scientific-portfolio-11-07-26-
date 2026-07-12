#!/usr/bin/env python3
"""Split sample assemblies into individual, traceable FASTA files."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


@dataclass(frozen=True)
class Record:
    identifier: str
    description: str
    sequence: str


def fasta_records(handle: TextIO) -> Iterator[Record]:
    header: str | None = None
    sequence: list[str] = []
    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield make_record(header, sequence)
            header = line[1:].strip()
            sequence = []
        elif header is None:
            raise ValueError(f"Sequence before FASTA header at line {line_number}")
        else:
            sequence.append(line)
    if header is not None:
        yield make_record(header, sequence)


def make_record(header: str, sequence: list[str]) -> Record:
    identifier = header.split(maxsplit=1)[0]
    if not identifier:
        raise ValueError("Empty FASTA identifier")
    return Record(identifier, header, "".join(sequence))


def safe_identifier(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    if not cleaned:
        raise ValueError(f"Identifier cannot be converted to a safe filename: {value!r}")
    return cleaned


def find_assembly(sample_dir: Path, candidates: list[str]) -> Path | None:
    for candidate in candidates:
        path = sample_dir / candidate
        if path.is_file():
            return path
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assemblies", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--assembly-name",
        action="append",
        dest="assembly_names",
        help="Candidate filename; may be repeated",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.assemblies.is_dir():
        raise SystemExit(f"Assembly directory not found: {args.assemblies}")
    names = args.assembly_names or ["assembly.fasta", "consensus.fasta"]
    args.output.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output / "contig_manifest.tsv"
    rows: list[dict[str, object]] = []

    for sample_dir in sorted(path for path in args.assemblies.iterdir() if path.is_dir()):
        assembly = find_assembly(sample_dir, names)
        if assembly is None:
            print(f"[WARN] No assembly for {sample_dir.name}")
            continue
        sample_id = safe_identifier(sample_dir.name)
        sample_out = args.output / sample_id
        sample_out.mkdir(parents=True, exist_ok=True)

        with assembly.open(encoding="utf-8") as handle:
            for record in fasta_records(handle):
                contig_id = safe_identifier(record.identifier)
                destination = sample_out / f"{contig_id}.fasta"
                if destination.exists() and not args.force:
                    raise FileExistsError(f"Refusing to overwrite {destination}; use --force")
                with destination.open("w", encoding="utf-8") as output:
                    output.write(f">{record.description}\n")
                    for start in range(0, len(record.sequence), 80):
                        output.write(record.sequence[start : start + 80] + "\n")
                rows.append(
                    {
                        "sample_id": sample_id,
                        "contig_id": record.identifier,
                        "length_bp": len(record.sequence),
                        "source_assembly": str(assembly.resolve()),
                        "contig_fasta": str(destination.resolve()),
                    }
                )

    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "contig_id", "length_bp", "source_assembly", "contig_fasta"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} contigs and {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
