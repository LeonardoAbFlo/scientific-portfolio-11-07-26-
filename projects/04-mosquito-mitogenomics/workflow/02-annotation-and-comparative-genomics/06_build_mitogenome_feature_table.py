#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path


GENUS_ABBR = {
    "Anopheles": "An.",
    "Culex": "Cx.",
    "Psorophora": "Ps.",
}

HEADER = [
    "Sample_ID",
    "Contig size (bp)",
    "PCGS (bp)",
    "AT%",
    "GC%",
    "A%",
    "T%",
    "G%",
    "C%",
    "AT SKEW",
    "GC SKEW",
]


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    workspace_root = script_dir.parent
    parser = argparse.ArgumentParser(
        description="Generate the mitochondrial genome features table for rrnS-corrected assemblies."
    )
    parser.add_argument(
        "--root-dir",
        type=Path,
        default=workspace_root / "results" / "annotation" / "corrections" / "medaka_mitofinder_rrnS",
        help="Directory that contains the genus subdirectories with corrected rrnS MitoFinder results.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=workspace_root / "results" / "reports" / "module2" / "mitogenome_features_table.tsv",
        help="Output TSV path.",
    )
    return parser.parse_args()


def format_percent(value: float) -> str:
    return f"{round(value, 1):.2f}".replace(".", ",")


def format_skew(value: float) -> str:
    return f"{value:.3f}".replace(".", ",")


def build_sample_id(sample_dir_name: str) -> str:
    barcode, taxon = sample_dir_name.split("-", 1)
    parts = taxon.split("_")
    genus = parts[0]
    species = parts[-1]
    if species.lower() in {"sp", "sp."}:
        species = "sp."
    return f"{GENUS_ABBR.get(genus, genus)} {species} ({barcode})"


def find_result_dir(sample_dir: Path) -> Path:
    internal_dirs = [path for path in sample_dir.iterdir() if path.is_dir() and path.name.endswith("_rrnS")]
    if len(internal_dirs) != 1:
        raise ValueError(f"Expected one internal _rrnS directory in {sample_dir}")
    result_dir = internal_dirs[0] / f"{internal_dirs[0].name}_MitoFinder_mitfi_Final_Results"
    if not result_dir.is_dir():
        raise ValueError(f"Missing MitoFinder results directory: {result_dir}")
    return result_dir


def get_unique_file(result_dir: Path, suffix: str) -> Path:
    matches = list(result_dir.glob(f"*{suffix}"))
    if len(matches) != 1:
        raise ValueError(f"Expected one file matching *{suffix} in {result_dir}")
    return matches[0]


def read_sequence(fasta_path: Path) -> str:
    return "".join(
        line.strip().upper()
        for line in fasta_path.read_text().splitlines()
        if line and not line.startswith(">")
    )


def calculate_pcgs_bp(gff_path: Path) -> int:
    total = 0
    for line in gff_path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        if "ATP8" in fields[8]:
            continue
        total += int(fields[4]) - int(fields[3]) + 1
    return total


def extract_barcode(sample_id: str) -> str:
    match = re.search(r"\(([^)]+)\)$", sample_id)
    if not match:
        raise ValueError(f"Could not extract barcode from {sample_id}")
    return match.group(1)


def build_rows(root_dir: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for genus_dir in sorted(path for path in root_dir.iterdir() if path.is_dir()):
        for sample_dir in sorted(path for path in genus_dir.iterdir() if path.is_dir()):
            result_dir = find_result_dir(sample_dir)
            fasta_path = get_unique_file(result_dir, "_mtDNA_contig.fasta")
            gff_path = get_unique_file(result_dir, "_mtDNA_contig.gff")

            sequence = read_sequence(fasta_path)
            counts = Counter(base for base in sequence if base in "ATGC")
            a_count, t_count, g_count, c_count = [counts[base] for base in "ATGC"]
            contig_size = a_count + t_count + g_count + c_count
            at_count = a_count + t_count
            gc_count = g_count + c_count

            rows.append(
                [
                    build_sample_id(sample_dir.name),
                    str(contig_size),
                    str(calculate_pcgs_bp(gff_path)),
                    format_percent(100 * at_count / contig_size),
                    format_percent(100 * gc_count / contig_size),
                    format_percent(100 * a_count / contig_size),
                    format_percent(100 * t_count / contig_size),
                    format_percent(100 * g_count / contig_size),
                    format_percent(100 * c_count / contig_size),
                    format_skew((a_count - t_count) / at_count if at_count else 0.0),
                    format_skew((g_count - c_count) / gc_count if gc_count else 0.0),
                ]
            )

    rows.sort(key=lambda row: extract_barcode(row[0]))
    return rows


def main() -> None:
    args = parse_args()
    root_dir = args.root_dir.resolve()
    output_path = args.output.resolve()
    rows = build_rows(root_dir)
    output_path.write_text(
        "\t".join(HEADER) + "\n" + "\n".join("\t".join(row) for row in rows) + "\n"
    )
    print(output_path)
    print(len(rows))


if __name__ == "__main__":
    main()
