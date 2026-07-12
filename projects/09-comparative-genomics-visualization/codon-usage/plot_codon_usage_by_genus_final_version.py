#!/usr/bin/env python3
"""Calculate RSCU tables and plots from mitochondrial GenBank files.

Input:
  <root>/<data-subdir>/<genus>/<sample_dir>/**/<gb-pattern>

Output:
  <root>/output/RSCU_plot/by_genus/<genus>/<sample_dir>/<barcode>/

The older ``mitofinder_trnI/barcode*`` layout is also accepted.
"""

from __future__ import annotations

from pathlib import Path
from collections import Counter, defaultdict
import re
import sys
import argparse

from Bio import SeqIO
from Bio.Data import CodonTable

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# Genetic code

GENETIC_CODE = "Invertebrate Mitochondrial"

aa3 = dict(
    A="Ala", R="Arg", N="Asn", D="Asp", C="Cys", Q="Gln", E="Glu",
    G="Gly", H="His", I="Ile", L="Leu", K="Lys", M="Met", F="Phe",
    P="Pro", S="Ser", T="Thr", W="Trp", Y="Tyr", V="Val", TER="TER"
)

tab = CodonTable.unambiguous_dna_by_name[GENETIC_CODE]
codon2aa = {c: aa3[a] for c, a in tab.forward_table.items()}
for stop in tab.stop_codons:
    codon2aa[stop] = "TER"

# Colors by third nucleotide
COLOR = dict(U="#b65ad0", T="#b65ad0", C="#1e90ff", A="#23a455", G="#e63946")


def codon_color(codon: str) -> str:
    """Return color based on the third nucleotide."""
    return COLOR.get(codon[2].replace("T", "U"), "#7f7f7f")


# Leucine and serine codon groups
leu1 = {"TTA", "TTG"}
leu2 = {"CTA", "CTC", "CTG", "CTT"}
ser1 = {"TCT", "TCC", "TCA", "TCG"}
ser2 = {"AGT", "AGC"}


def fam_label(codon: str, aa: str) -> str:
    """Separate leucine and serine codons into their two groups."""
    if aa == "Leu":
        return "Leu1" if codon in leu1 else "Leu2"
    if aa == "Ser":
        return "Ser1" if codon in ser1 else "Ser2"
    return aa


AA_ORDER = [
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly",
    "His", "Ile", "Leu1", "Leu2", "Lys", "Met", "Phe", "Pro",
    "Ser1", "Ser2", "Thr", "Trp", "Tyr", "Val"
]


# Path helpers

def extract_barcode_name(text: str) -> str:
    """Extract ``barcodeNN`` or return the filename stem."""
    m = re.search(r"(barcode\d+)", text)
    return m.group(1) if m else Path(text).stem


def genus_from_path(file_path: Path, data_dir: Path) -> str:
    """
    <data>/<genus>/<sample_dir>/.../<file> -> genus
    """
    try:
        rel = file_path.resolve().relative_to(data_dir.resolve())
        if len(rel.parts) >= 1:
            return rel.parts[0]
    except Exception:
        pass
    return "unknown"


def sample_name_from_path(file_path: Path, data_dir: Path) -> str:
    """
    <data>/<genus>/<sample_dir>/.../<file> -> sample_dir
    """
    try:
        rel = file_path.resolve().relative_to(data_dir.resolve())
        if len(rel.parts) >= 2:
            return rel.parts[1]
    except Exception:
        pass
    return file_path.parent.name


def parse_sample_dir(sample_dir_name: str) -> tuple[str, str]:
    """
    'F3_92-Anopheles_nyssorhynchus_triannulatus'
      -> ('F3_92', 'Anopheles nyssorhynchus triannulatus')
    If '-' is absent, use the directory name as the species label.
    """
    if "-" in sample_dir_name:
        sample_id, mosq = sample_dir_name.split("-", 1)
    else:
        sample_id, mosq = "", sample_dir_name

    mosq = mosq.replace("_", " ")
    mosq = re.sub(r"\s+", " ", mosq).strip()
    sample_id = sample_id.strip()

    return sample_id, mosq


def math_italic(text: str) -> str:
    """Return italic text formatted for Matplotlib mathtext."""
    # Escape characters used by mathtext
    esc = text.replace("\\", r"\\").replace("{", r"\{").replace("}", r"\}")
    # Preserve spaces in mathtext
    esc = esc.replace(" ", r"\ ")
    return rf"$\it{{{esc}}}$"


def build_plot_title(mosq_name: str, sample_id: str) -> str:
    """Build a title with italic species name and sample ID."""
    mosq_it = math_italic(mosq_name)
    return f"{mosq_it} ({sample_id})" if sample_id else mosq_it


# RSCU calculation

def rscu_from_genbank(
    gb_path: Path,
    barcode: str,
    out_prefix: Path,
    plot_title: str
) -> tuple[Path, Path]:
    """Calculate RSCU and create a plot for one mitochondrial GenBank file."""
    gb_path = Path(gb_path)
    out_prefix = Path(out_prefix)

    # Count codons in CDS features
    cnt = Counter()
    any_record = False
    for record in SeqIO.parse(str(gb_path), "genbank"):
        any_record = True
        for feat in record.features:
            if feat.type == "CDS":
                seq = feat.extract(record.seq).upper()
                # Use complete codons only
                for i in range(0, len(seq) - 2, 3):
                    cod = str(seq[i:i + 3])
                    if cod in codon2aa:
                        cnt[cod] += 1

    if not any_record:
        raise ValueError(f"GenBank file is empty or cannot be parsed: {gb_path}")

    # Calculate RSCU and include stop codons in the table
    aa2cod = defaultdict(list)
    for cod, aa in codon2aa.items():
        aa2cod[aa].append(cod)

    rows = []
    for aa in sorted(aa2cod):
        cods = aa2cod[aa]
        total = sum(cnt[c] for c in cods)
        n = len(cods)
        for cod in cods:
            obs = cnt[cod]
            rscu = (obs * n) / total if total else 0.0
            rows.append([aa, cod, obs, round(float(rscu), 2)])

    df = pd.DataFrame(rows, columns=["AminoAcid", "Codon", "Count", "RSCU"])
    df["Group"] = [fam_label(c, a) for c, a in zip(df["Codon"], df["AminoAcid"])]

    # Exclude stop codons from the plot
    df_plot = df[df["Group"] != "TER"].copy()

    # Codons grouped by amino acid
    grp2cod = defaultdict(list)
    for cod, grp in zip(df_plot["Codon"], df_plot["Group"]):
        grp2cod[grp].append(cod)

    # Plot
    plt.rcParams.update({"font.size": 11, "font.family": "sans-serif"})
    fig = plt.figure(figsize=(16, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1.4], hspace=0.05)

    # Stacked bars by amino acid group
    ax = fig.add_subplot(gs[0])
    x = np.arange(len(AA_ORDER))
    bottom = np.zeros_like(x, dtype=float)
    bar_w = 0.8

    for cod in sorted(df_plot["Codon"].unique()):
        grp = df_plot.loc[df_plot["Codon"] == cod, "Group"].iloc[0]
        if grp not in AA_ORDER:
            continue
        pos = AA_ORDER.index(grp)
        height = float(df_plot.loc[df_plot["Codon"] == cod, "RSCU"].iloc[0])

        ax.bar(
            pos,
            height,
            bar_w,
            bottom=bottom[pos],
            color=codon_color(cod),
            edgecolor="white",
            linewidth=0.4,
        )
        bottom[pos] += height

    ax.set_xticks(x)
    ax.set_xticklabels(AA_ORDER, rotation=0, ha="center")
    ax.tick_params(axis="x", pad=6)
    ax.set_ylabel("RSCU")

    # Species and sample title
    ax.set_title(plot_title)

    # Codon labels below the bars
    cod_ax = fig.add_subplot(gs[1], sharex=ax)
    cod_ax.set_xlim(-0.5, len(AA_ORDER) - 0.5)
    cod_ax.set_ylim(0, 1)
    cod_ax.axis("off")

    box_h, v_gap, y0 = 0.18, 0.02, 0.06
    for j, grp in enumerate(AA_ORDER):
        codons = sorted(grp2cod.get(grp, []))[::-1]  # de arriba abajo
        for i, cod in enumerate(codons):
            y = y0 + i * (box_h + v_gap)
            color = codon_color(cod)
            cod_ax.add_patch(
                Rectangle(
                    (j - bar_w / 2, y),
                    bar_w,
                    box_h,
                    facecolor=color,
                    edgecolor="k",
                    linewidth=0,
                )
            )
            cod_ax.text(
                j,
                y + box_h / 2,
                cod,
                ha="center",
                va="center",
                fontsize=8,
                color="white" if color != "#1e90ff" else "black",
            )

    # Save table and plot
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig_path = out_prefix.parent / f"{barcode}_rscu.png"
    csv_path = out_prefix.parent / f"{barcode}_rscu_table.csv"

    fig.tight_layout()
    fig.savefig(fig_path, dpi=600, bbox_inches="tight")
    plt.close(fig)

    df.to_csv(csv_path, index=False)

    print(f"[OK] {barcode}  →  {fig_path} , {csv_path}")
    return fig_path, csv_path


# Input discovery

def find_genbank_in_results_layout(data_dir: Path, gb_pattern: str) -> list[tuple[str, str, str, Path]]:
    """Find GenBank files grouped by genus and sample directory."""
    pairs: list[tuple[str, str, str, Path]] = []
    for gb_path in sorted(Path(data_dir).rglob(gb_pattern)):
        gb_path = gb_path.resolve()
        genus = genus_from_path(gb_path, data_dir)
        sample = sample_name_from_path(gb_path, data_dir)

        # Read the barcode from the file or its parent directory
        barcode = extract_barcode_name(gb_path.name)
        if barcode == Path(gb_path.name).stem:
            barcode = extract_barcode_name(sample)
        if barcode == Path(sample).stem:
            barcode = extract_barcode_name(gb_path.parent.name)

        pairs.append((genus, sample, barcode, gb_path))
    return pairs


def find_genbank_in_old_mitof_layout(mitof_root: Path, gb_pattern: str) -> list[tuple[str, str, str, Path]]:
    """Find GenBank files in the older ``mitofinder_trnI`` layout."""
    pairs: list[tuple[str, str, str, Path]] = []
    mitof_root = Path(mitof_root)
    for bc_dir in sorted(mitof_root.glob("barcode*_trnI")):
        barcode = extract_barcode_name(bc_dir.name)  # barcodeNN
        gb_candidates = list(bc_dir.rglob(gb_pattern))
        if not gb_candidates:
            print(f"[WARN] No files matching {gb_pattern} in {bc_dir}", file=sys.stderr)
            continue
        gb_path = gb_candidates[0].resolve()
        pairs.append(("unknown", bc_dir.name, barcode, gb_path))
    return pairs


# Main

def main():
    p = argparse.ArgumentParser(
        description="Calculate RSCU by genus and sample from a results directory."
    )
    p.add_argument("--root", default="/path/to/mosquitos/results_p1", help="Project root")
    p.add_argument("--data-subdir", default="CORRECTIONS_medaka_mitofinder_rrnS", help="Subdirectory containing <genus>/<sample_dir>")
    p.add_argument("--outdir", default=None, help="Output directory (default: <root>/output/RSCU_plot)")
    p.add_argument("--gb-pattern", default="*_mtDNA_contig.gb", help="GenBank filename pattern")
    p.add_argument("--mitof-subdir", default="mitofinder_trnI", help="Older layout under <root>/<mitof-subdir>")
    p.add_argument("--species-name", default=None, help="Override the species name; sample ID is read from the directory")
    args = p.parse_args()

    root = Path(args.root).expanduser().resolve()
    data_dir = (root / args.data_subdir).resolve()

    output_dir = Path(args.outdir).expanduser().resolve() if args.outdir else (root / "output" / "RSCU_plot").resolve()
    out_base_by_genus = output_dir / "by_genus"
    out_base_by_barcode = output_dir / "by_barcode"

    pairs: list[tuple[str, str, str, Path]] = []
    mode = None

    if data_dir.is_dir():
        mode = "results_layout"
        pairs = find_genbank_in_results_layout(data_dir, args.gb_pattern)
    else:
        mode = "old_mitof_layout"
        mitof_root = (root / args.mitof_subdir).resolve()
        if mitof_root.is_dir():
            pairs = find_genbank_in_old_mitof_layout(mitof_root, args.gb_pattern)
        else:
            raise SystemExit(f"[ERROR] Neither {data_dir} nor {mitof_root} exists")

    if not pairs:
        raise SystemExit(f"[ERROR] No GenBank files matching {args.gb_pattern} in {mode} mode")

    print(f"[INFO] Mode: {mode} | Found {len(pairs)} GenBank files.")
    print(f"[INFO] Output directory: {output_dir}")

    for genus, sample_dir, barcode, gb_path in pairs:
        # Read sample ID and species from the directory name
        sample_id, mosq_from_sample = parse_sample_dir(sample_dir)

        # Use the command-line species name when provided
        mosq_name = args.species_name.strip() if args.species_name else mosq_from_sample

        # Species and sample title
        plot_title = build_plot_title(mosq_name, sample_id)

        try:
            if mode == "results_layout":
                # .../output/RSCU_plot/by_genus/<genus>/<sample>/<barcode>/
                out_prefix = out_base_by_genus / genus / sample_dir / barcode / barcode
            else:
                # .../output/RSCU_plot/by_barcode/<barcode>/
                out_prefix = out_base_by_barcode / barcode / barcode

            print(f"[INFO] Processing {barcode} ({gb_path})")
            rscu_from_genbank(gb_path, barcode=barcode, out_prefix=out_prefix, plot_title=plot_title)

        except Exception as e:
            print(f"[ERROR] {barcode}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
