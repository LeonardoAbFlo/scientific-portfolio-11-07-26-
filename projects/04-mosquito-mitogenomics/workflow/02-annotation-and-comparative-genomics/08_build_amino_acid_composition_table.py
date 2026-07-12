#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


GENUS_ABBR = {
    "Anopheles": "An.",
    "Culex": "Cx.",
    "Psorophora": "Ps.",
}

GENE_ORDER = [
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6",
]

PCG_GENES = set(GENE_ORDER)

AA_SYMBOL_TO_NAME = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
}

AA_ORDER = [
    "Ala",
    "Arg",
    "Asn",
    "Asp",
    "Cys",
    "Gln",
    "Glu",
    "Gly",
    "His",
    "Ile",
    "Leu",
    "Lys",
    "Met",
    "Phe",
    "Pro",
    "Ser",
    "Thr",
    "Trp",
    "Tyr",
    "Val",
]

PREFERRED_SAMPLE_ORDER = [
    "F3_92",
    "F3_95",
    "F4_90",
    "F1_83",
    "F3_93",
    "F4_93",
    "F1_81",
    "F1_84",
    "F3_90",
    "F3_89",
]


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    workspace_root = script_dir.parent
    parser = argparse.ArgumentParser(
        description="Generate Table S4 amino acid composition percentages from rrnS-corrected MitoFinder results."
    )
    parser.add_argument(
        "--root-dir",
        type=Path,
        default=workspace_root / "results" / "annotation" / "corrections" / "medaka_mitofinder_rrnS",
        help="Directory that contains the corrected rrnS MitoFinder result folders.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=workspace_root / "results" / "reports" / "module2" / "table_s4_amino_acid_composition.tsv",
        help="Output TSV path.",
    )
    return parser.parse_args()


def find_result_dir(sample_dir: Path) -> Path:
    internal_dirs = [path for path in sample_dir.iterdir() if path.is_dir() and path.name.endswith("_rrnS")]
    if len(internal_dirs) != 1:
        raise ValueError(f"Expected one internal _rrnS directory in {sample_dir}")
    result_dir = internal_dirs[0] / f"{internal_dirs[0].name}_MitoFinder_mitfi_Final_Results"
    if not result_dir.is_dir():
        raise ValueError(f"Missing MitoFinder results directory: {result_dir}")
    return result_dir


def get_aa_fasta(result_dir: Path) -> Path:
    matches = list(result_dir.glob("*_mtDNA_contig_genes_AA.fasta"))
    if len(matches) != 1:
        raise ValueError(f"Expected one *_mtDNA_contig_genes_AA.fasta file in {result_dir}")
    return matches[0]


def parse_fasta(fasta_path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                sequences[header] = "".join(chunks).upper()
            header = line[1:]
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        sequences[header] = "".join(chunks).upper()
    return sequences


def get_gene_name(header: str) -> str:
    return header.split("@")[-1]


def build_species_label(sample_dir_name: str) -> str:
    barcode, taxon = sample_dir_name.split("-", 1)
    parts = taxon.split("_")
    genus = parts[0]
    species = parts[-1]
    if species.lower() in {"sp", "sp."}:
        return f"{genus} sp. ({barcode})"
    return f"{GENUS_ABBR.get(genus, genus)} {species} ({barcode})"


def barcode_sort_key(species_label: str) -> tuple[int, str]:
    barcode = species_label.rsplit("(", 1)[-1].rstrip(")")
    if barcode in PREFERRED_SAMPLE_ORDER:
        return (PREFERRED_SAMPLE_ORDER.index(barcode), barcode)
    return (len(PREFERRED_SAMPLE_ORDER), barcode)


def count_amino_acids(aa_fasta_path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    sequences = parse_fasta(aa_fasta_path)
    for header, sequence in sequences.items():
        gene = get_gene_name(header)
        if gene not in PCG_GENES:
            continue
        for aa_symbol in sequence:
            if aa_symbol == "*":
                continue
            amino_acid = AA_SYMBOL_TO_NAME.get(aa_symbol)
            if amino_acid is None:
                raise ValueError(f"Unexpected amino acid symbol '{aa_symbol}' in {aa_fasta_path}")
            counts[amino_acid] += 1
    return counts


def format_percentage(value: float) -> str:
    return f"{value:.3f}"


def build_rows(root_dir: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for genus_dir in sorted(path for path in root_dir.iterdir() if path.is_dir()):
        for sample_dir in sorted(path for path in genus_dir.iterdir() if path.is_dir()):
            result_dir = find_result_dir(sample_dir)
            aa_fasta = get_aa_fasta(result_dir)
            counts = count_amino_acids(aa_fasta)
            total = sum(counts.values())
            if total == 0:
                raise ValueError(f"No amino acid residues counted in {aa_fasta}")

            row = {"Species": build_species_label(sample_dir.name)}
            for amino_acid in AA_ORDER:
                row[amino_acid] = format_percentage(counts[amino_acid] * 100 / total)
            rows.append(row)
    rows.sort(key=lambda row: barcode_sort_key(row["Species"]))
    return rows


def write_tsv(output_path: Path, rows: list[dict[str, str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Species", *AA_ORDER], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    rows = build_rows(args.root_dir.resolve())
    write_tsv(args.output.resolve(), rows)
    print(args.output.resolve())
    print(len(rows))


if __name__ == "__main__":
    main()
