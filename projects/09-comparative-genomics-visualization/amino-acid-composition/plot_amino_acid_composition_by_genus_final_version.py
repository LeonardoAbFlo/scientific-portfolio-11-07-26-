#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import argparse
import re
from collections import Counter
from typing import List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from Bio import SeqIO

# Plot settings
FIGSIZE = (13, 5)
DPI = 600

# Line settings
LINEWIDTH = 1.2
ALPHA_LINE = 0.40

# Point settings
MARKERSIZE = 4
ALPHA_MARKER = 0.6

TITLE_SIZE = 14
LABEL_SIZE = 12
TICK_SIZE = 10
LEGEND_SIZE = 8

# Colors
BASE_HEX = [
    "#ffbd00",
    "#ff5400",
    "#ff0054",
    "#9e0059",
    "#390099",
    "#42bfdd",
    "#06d6a0",
]

def make_extended_palette(n: int, base_hex=BASE_HEX):
    """
    Create n colors by smoothly interpolating the base palette.
    Returns a list of RGBA tuples usable directly in matplotlib.
    """
    if n <= 1:
        return [mpl.colors.to_rgba(base_hex[0])]

    cmap = mpl.colors.LinearSegmentedColormap.from_list("custom_extended", base_hex)
    vals = np.linspace(0.0, 1.0, n)
    return [cmap(v) for v in vals]

# Amino acids
AA1_TO_AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}

AA_ORDER_1L = list(AA1_TO_AA3.keys())
AA_ORDER_3L = [AA1_TO_AA3[a] for a in AA_ORDER_1L]

# Utilities
def parse_sample_dir(sample_dir: str) -> tuple[str, str]:
    """
    sample_dir:
      'F1_83-Culex_quinquefasciatus'
    returns:
      sample_id = 'F1_83'
      species   = 'Culex quinquefasciatus'
    """
    if "-" in sample_dir:
        left, right = sample_dir.split("-", 1)
        sample_id = left.strip()
        species = right.strip()
    else:
        sample_id = sample_dir.strip()
        species = sample_dir.strip()

    species = species.replace("_", " ")
    species = re.sub(r"\s+", " ", species).strip()
    return sample_id, species


def italicize_species_mathtext(species: str) -> str:
    """
    Italic species name using matplotlib mathtext.
    E.g. 'Culex quinquefasciatus' -> r'$\\it{Culex\\ quinquefasciatus}$'
    """
    species_mt = species.replace(" ", r"\ ")
    return rf"$\it{{{species_mt}}}$"


def make_label(sample_dir: str) -> str:
    """
    'F1_83-Culex_quinquefasciatus' -> '$\\it{Culex\\ quinquefasciatus}$ (F1_83)'
    """
    sample_id, species = parse_sample_dir(sample_dir)
    species_it = italicize_species_mathtext(species)
    return f"{species_it} ({sample_id})"


def load_aa_sequences(fasta: Path) -> List[str]:
    return [str(rec.seq).upper() for rec in SeqIO.parse(str(fasta), "fasta")]


def aa_counts(seqs: List[str]) -> Tuple[Counter, int]:
    cnt = Counter()
    total = 0
    for s in seqs:
        for aa in s:
            if aa in AA1_TO_AA3:
                cnt[aa] += 1
                total += 1
    return cnt, total


def counts_to_series(cnt: Counter, total: int) -> pd.Series:
    return pd.Series(
        {AA1_TO_AA3[a]: (cnt.get(a, 0) / total * 100 if total else 0.0)
         for a in AA_ORDER_1L}
    ).reindex(AA_ORDER_3L)

# Plot
def plot_aa_composition(df: pd.DataFrame, title: str, out_png: Path):
    fig, ax = plt.subplots(figsize=FIGSIZE)

    # Background and grid
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_axisbelow(True)

    # Keep labels in the same order
    labels = sorted(df["Label"].unique())
    colors = make_extended_palette(len(labels))

    for i, label in enumerate(labels):
        sub = df[df["Label"] == label]

        # Draw the connecting line
        ax.plot(
            sub["AA"], sub["Percent"],
            label=label,
            linewidth=LINEWIDTH,
            alpha=ALPHA_LINE,
            color=colors[i],
            zorder=2
        )

        # Draw the points
        ax.plot(
            sub["AA"], sub["Percent"],
            linestyle="None",
            marker="o",
            markersize=MARKERSIZE,
            alpha=ALPHA_MARKER,
            color=colors[i],
            zorder=3
        )

    ax.set_title(title, fontsize=TITLE_SIZE)
    ax.set_xlabel("Amino acids", fontsize=LABEL_SIZE)
    ax.set_ylabel("Percentage (%)", fontsize=LABEL_SIZE)
    ax.tick_params(labelsize=TICK_SIZE)

    ax.set_ylim(bottom=0)
    ax.grid(True, which="major", linestyle="-", linewidth=0.8, color="#E6E6E6")

    # Legend
    ax.legend(
        frameon=False,
        fontsize=LEGEND_SIZE,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0,
        handlelength=1.4,
        handletextpad=0.6,
        labelspacing=0.35,
        markerscale=0.9
    )

    # Leave space for the legend
    fig.tight_layout(rect=[0, 0, 0.80, 1])

    # Save the plot
    fig.savefig(out_png, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_png}")

# Main
def main():
    ap = argparse.ArgumentParser(
        description="Amino acid composition per genus and global (species label includes sample ID)."
    )
    ap.add_argument("--root", default="/path/to/mosquitos/results_p1")
    ap.add_argument("--data-subdir", default="CORRECTIONS_medaka_mitofinder_rrnS")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    root = Path(args.root).resolve()
    data_dir = root / args.data_subdir
    if not data_dir.is_dir():
        raise SystemExit(f"[ERROR] Not found: {data_dir}")

    outdir = Path(args.outdir).resolve() if args.outdir else root / "output" / "aa_composition"
    outdir.mkdir(parents=True, exist_ok=True)
    by_genus_dir = outdir / "by_genus"
    by_genus_dir.mkdir(exist_ok=True)

    rows_all = []

    # By genus
    for genus_dir in sorted(p for p in data_dir.iterdir() if p.is_dir()):
        genus = genus_dir.name
        rows_genus = []

        for sample_dir in sorted(p for p in genus_dir.iterdir() if p.is_dir()):
            fasta_files = list(sample_dir.rglob("*_final_genes_AA.fasta"))
            if not fasta_files:
                continue

            label = make_label(sample_dir.name)
            seqs = load_aa_sequences(fasta_files[0])
            cnt, total = aa_counts(seqs)
            s = counts_to_series(cnt, total)

            for aa, pct in s.items():
                row = {"Genus": genus, "Label": label, "AA": aa, "Percent": pct}
                rows_genus.append(row)
                rows_all.append(row)

        if not rows_genus:
            continue

        df_genus = pd.DataFrame(rows_genus)
        df_genus["AA"] = pd.Categorical(df_genus["AA"], AA_ORDER_3L, ordered=True)
        df_genus = df_genus.sort_values(["AA", "Label"])

        df_genus.to_csv(outdir / f"{genus}_aa_composition.csv", index=False)

        plot_aa_composition(
            df_genus,
            title=f"Amino acid composition – {genus.capitalize()}",
            out_png=by_genus_dir / f"{genus}_aa_composition.png"
        )

    # Global
    if rows_all:
        df_all = pd.DataFrame(rows_all)
        df_all["AA"] = pd.Categorical(df_all["AA"], AA_ORDER_3L, ordered=True)
        df_all = df_all.sort_values(["AA", "Label"])

        df_all.to_csv(outdir / "all_aa_composition.csv", index=False)

        plot_aa_composition(
            df_all,
            title="Amino acid composition – Culicidae mitogenomes",
            out_png=outdir / "all_aa_composition.png"
        )
    else:
        print("[WARN] No rows generated (missing *_final_genes_AA.fasta?)")


if __name__ == "__main__":
    main()
