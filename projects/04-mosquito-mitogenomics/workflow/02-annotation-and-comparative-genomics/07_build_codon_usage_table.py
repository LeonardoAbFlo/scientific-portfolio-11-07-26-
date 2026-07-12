#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
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
    "TER",
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
    "*": "TER",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
}

STOP_CODONS = {"TAA", "TAG"}


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    workspace_root = script_dir.parent
    parser = argparse.ArgumentParser(
        description="Generate Table S3 codon usage summaries for rrnS-corrected mitochondrial genomes."
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
        default=workspace_root / "results" / "reports" / "module2" / "table_s3_codon_usage.tsv",
        help="Output TSV path.",
    )
    parser.add_argument(
        "--long-output",
        type=Path,
        default=workspace_root / "results" / "reports" / "module2" / "table_s3_codon_usage_long.tsv",
        help="Output TSV path for the long by-gene codon table.",
    )
    parser.add_argument(
        "--min-fraction",
        type=float,
        default=0.5,
        help="Minimum fraction of the top codon count required to keep additional codons in a cell.",
    )
    parser.add_argument(
        "--max-codons-per-cell",
        type=int,
        default=3,
        help="Maximum number of codons to report per amino acid and sample.",
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


def get_gene_fastas(result_dir: Path) -> tuple[Path, Path]:
    nt_matches = list(result_dir.glob("*_mtDNA_contig_genes_NT.fasta"))
    aa_matches = list(result_dir.glob("*_mtDNA_contig_genes_AA.fasta"))
    if len(nt_matches) != 1:
        raise ValueError(f"Expected one *_mtDNA_contig_genes_NT.fasta file in {result_dir}")
    if len(aa_matches) != 1:
        raise ValueError(f"Expected one *_mtDNA_contig_genes_AA.fasta file in {result_dir}")
    return nt_matches[0], aa_matches[0]


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


def build_species_label(sample_dir_name: str) -> str:
    barcode, taxon = sample_dir_name.split("-", 1)
    parts = taxon.split("_")
    genus = parts[0]
    species = parts[-1]
    if species.lower() in {"sp", "sp."}:
        return f"{GENUS_ABBR.get(genus, genus)} sp. ({barcode})"
    return f"{GENUS_ABBR.get(genus, genus)} {species} ({barcode})"


def barcode_sort_key(species_label: str) -> tuple[int, str]:
    barcode = species_label.rsplit("(", 1)[-1].rstrip(")")
    if barcode in PREFERRED_SAMPLE_ORDER:
        return (PREFERRED_SAMPLE_ORDER.index(barcode), barcode)
    return (len(PREFERRED_SAMPLE_ORDER), barcode)


def gene_sort_key(gene: str) -> tuple[int, str]:
    if gene in GENE_ORDER:
        return (GENE_ORDER.index(gene), gene)
    return (len(GENE_ORDER), gene)


def get_gene_name(header: str) -> str:
    return header.split("@")[-1]


def normalize_aa_sequence(codons: list[str], amino_acids: str) -> str:
    if len(codons) == len(amino_acids):
        return amino_acids
    if len(codons) == len(amino_acids) + 1 and codons[-1] in STOP_CODONS:
        return amino_acids + "*"
    raise ValueError(
        "Gene NT/AA length mismatch: "
        f"{len(codons)} codons versus {len(amino_acids)} amino acids"
    )


def collect_codon_counts(
    nt_fasta_path: Path,
    aa_fasta_path: Path,
) -> tuple[dict[str, Counter[str]], dict[str, dict[str, Counter[str]]]]:
    summary_counts: dict[str, Counter[str]] = defaultdict(Counter)
    gene_counts: dict[str, dict[str, Counter[str]]] = defaultdict(lambda: defaultdict(Counter))
    nt_sequences = {
        get_gene_name(header): sequence
        for header, sequence in parse_fasta(nt_fasta_path).items()
    }
    aa_sequences = {
        get_gene_name(header): sequence
        for header, sequence in parse_fasta(aa_fasta_path).items()
    }

    for gene, nt_sequence in nt_sequences.items():
        if gene not in PCG_GENES:
            continue
        if gene not in aa_sequences:
            raise ValueError(f"Missing AA sequence for gene {gene} in {aa_fasta_path}")

        codons = [nt_sequence[start:start + 3] for start in range(0, len(nt_sequence) - 2, 3)]
        amino_acids = normalize_aa_sequence(codons, aa_sequences[gene])

        for codon, aa_symbol in zip(codons, amino_acids):
            amino_acid = AA_SYMBOL_TO_NAME.get(aa_symbol)
            if amino_acid is None:
                continue
            summary_counts[amino_acid][codon] += 1
            gene_counts[gene][amino_acid][codon] += 1
    return summary_counts, gene_counts


def select_codons(
    codon_counts: Counter[str],
    min_fraction: float,
    max_codons_per_cell: int,
) -> list[tuple[str, int]]:
    if not codon_counts:
        return []
    ranked = sorted(codon_counts.items(), key=lambda item: (-item[1], item[0]))
    top_count = ranked[0][1]
    selected = [
        (codon, count)
        for codon, count in ranked
        if count >= top_count * min_fraction
    ]
    return selected[:max_codons_per_cell]


def summarize_codons(codon_counts: Counter[str], min_fraction: float, max_codons_per_cell: int) -> str:
    selected = select_codons(codon_counts, min_fraction, max_codons_per_cell)
    return "; ".join(f"{codon} ({count})" for codon, count in selected)


def build_tables(
    root_dir: Path,
    min_fraction: float,
    max_codons_per_cell: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    summary_rows: list[dict[str, str]] = []
    long_rows: list[dict[str, str]] = []
    for genus_dir in sorted(path for path in root_dir.iterdir() if path.is_dir()):
        for sample_dir in sorted(path for path in genus_dir.iterdir() if path.is_dir()):
            result_dir = find_result_dir(sample_dir)
            nt_fasta, aa_fasta = get_gene_fastas(result_dir)
            sample_counts, gene_counts = collect_codon_counts(nt_fasta, aa_fasta)
            barcode = sample_dir.name.split("-", 1)[0]
            species_label = build_species_label(sample_dir.name)

            row = {"Species": species_label}
            for amino_acid in AA_ORDER:
                row[amino_acid] = summarize_codons(sample_counts.get(amino_acid, Counter()), min_fraction, max_codons_per_cell)
            summary_rows.append(row)

            selected_codons_by_aa = {
                amino_acid: {codon for codon, _ in select_codons(codon_counts, min_fraction, max_codons_per_cell)}
                for amino_acid, codon_counts in sample_counts.items()
            }
            for gene in sorted(gene_counts, key=gene_sort_key):
                for amino_acid in AA_ORDER:
                    codon_counts = gene_counts[gene].get(amino_acid)
                    if not codon_counts:
                        continue
                    for codon, count in sorted(codon_counts.items(), key=lambda item: (-item[1], item[0])):
                        long_rows.append(
                            {
                                "Species": species_label,
                                "Barcode": barcode,
                                "Gene": gene,
                                "AminoAcid": amino_acid,
                                "Codon": codon,
                                "Gene_Codon_Count": str(count),
                                "Sample_Codon_Total": str(sample_counts[amino_acid][codon]),
                                "Selected_For_Summary": "yes" if codon in selected_codons_by_aa.get(amino_acid, set()) else "no",
                            }
                        )
    summary_rows.sort(key=lambda row: barcode_sort_key(row["Species"]))
    long_rows.sort(
        key=lambda row: (
            barcode_sort_key(row["Species"]),
            gene_sort_key(row["Gene"]),
            AA_ORDER.index(row["AminoAcid"]),
            -int(row["Gene_Codon_Count"]),
            row["Codon"],
        )
    )
    return summary_rows, long_rows


def write_tsv(output_path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    summary_rows, long_rows = build_tables(
        args.root_dir.resolve(),
        args.min_fraction,
        args.max_codons_per_cell,
    )
    write_tsv(args.output.resolve(), summary_rows, ["Species", *AA_ORDER])
    write_tsv(
        args.long_output.resolve(),
        long_rows,
        [
            "Species",
            "Barcode",
            "Gene",
            "AminoAcid",
            "Codon",
            "Gene_Codon_Count",
            "Sample_Codon_Total",
            "Selected_For_Summary",
        ],
    )
    print(args.output.resolve())
    print(len(summary_rows))
    print(args.long_output.resolve())
    print(len(long_rows))


if __name__ == "__main__":
    main()
