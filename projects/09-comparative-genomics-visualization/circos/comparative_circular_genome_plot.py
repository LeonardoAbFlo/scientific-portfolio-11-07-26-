#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Create a comparative circular plot for mosquito mitogenomes.

BLASTN is used for nucleotide comparisons. The outer ring shows CDS in black,
tRNA in green, and rRNA in blue. Input samples are found under
``<root>/data/<genus>/<sample_dir>``. The reference sample can be configured from
the command line.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import tempfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge
from Bio import SeqIO
from Bio.Seq import Seq


# Samples and colors
SELECTED_SAMPLE_COLORS = OrderedDict([
    ("Anopheles triannulatus (F3_92)", "#ff9896"),
    ("Culex quinquefasciatus (F3_93)", "#a0522d"),
    ("Culex iolambis (F3_95)", "#9467bd"),
    ("Culex sp. (F4_93)", "#8c564b"),
    ("Psorophora ferox (F3_89)", "#2ca02c"),
    ("Psorophora sp. (F3_90)", "#98df8a"),
    ("Coquilletidia nigricans (F1_82)", "#1f77b4"),
    ("Coquillettidia venezuelensis (F1_85)", "#6baed6"),
    ("Coquilletidia sp. (F3_91)", "#d62728"),
    ("Ochlerotatus serratus (F3_96)", "#aec7e8"),
    ("Ochlerotatus fulvus (F3_94)", "#ffbb78"),
    ("Limatus flavisetosus (F1_87)", "#ff7f0e"),
    ("Mansonia humeralis (F5_47)", "#636363"),
    ("Mansonia sp. (F4_89)", "#0e8c9c"),
    ("Mansonia sp. (F5_46)", "#17becf"),
    ("Mansonia sp. (F5_48)", "#8fe3eb"),
])

REFERENCE_LABEL = "Anopheles triannulatus (F3_92)"
COMPARISON_LABELS = [x for x in SELECTED_SAMPLE_COLORS.keys() if x != REFERENCE_LABEL]


# Gene colors
FEATURE_COLORS = {
    "CDS": "#000000",    # PCGs / protein-coding genes
    "tRNA": "#2ca02c",   # transfer RNAs
    "rRNA": "#1f77b4",   # ribosomal RNAs
}

FEATURE_LEGEND_LABELS = {
    "CDS": "PCGs",
    "tRNA": "tRNAs",
    "rRNA": "rRNAs",
}

BACKGROUND_TRACK_COLOR = "#f0f0f0"
GENE_LABEL_COLOR = "#111111"
INNER_TICK_COLOR = "#707070"
CENTER_TEXT_COLOR = "#111111"
LEADER_LINE_COLOR = "#6a6a6a"


# Data structures
@dataclass
class Feature:
    start: int
    end: int
    strand: str
    ftype: str
    label: str


@dataclass
class GenomeRecord:
    display_label: str
    sample_code: str
    root: Path
    sample_dir: Path
    seq_path: Path
    seq_id: str
    seq: Seq
    gb_path: Optional[Path] = None
    gff_path: Optional[Path] = None
    features: Optional[list[Feature]] = None

    @property
    def length(self) -> int:
        return len(self.seq)


# File and sample helpers
def extract_sample_code(display_label: str) -> str:
    m = re.search(r"\((F\d+_\d+)\)", display_label)
    if not m:
        raise ValueError(f"Could not extract an F?_? sample code from: {display_label}")
    return m.group(1)


def normalize_text(s: str) -> str:
    s = s.lower().strip()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def preferred_tokens_from_label(display_label: str) -> list[str]:
    species_text = re.sub(r"\s*\(F\d+_\d+\)$", "", display_label).strip()
    return [t for t in normalize_text(species_text).split("_") if t not in {"sp"}]


def get_best_sample_dir(roots: list[Path], display_label: str) -> Path:
    code = extract_sample_code(display_label)
    desired_tokens = preferred_tokens_from_label(display_label)
    candidates: list[Path] = []

    for root in roots:
        data_dir = root / "data"
        if data_dir.is_dir():
            candidates.extend(sorted(data_dir.glob(f"*/{code}-*")))

    if not candidates:
        raise FileNotFoundError(
            f"Directory not found for {display_label}. Searched for {code}-* in: "
            + ", ".join(str(r / "data") for r in roots)
        )

    if len(candidates) == 1:
        return candidates[0]

    def score(path: Path) -> tuple[int, int]:
        name = normalize_text(path.name)
        hits = sum(tok in name for tok in desired_tokens)
        return hits, len(name)

    return sorted(candidates, key=score, reverse=True)[0]


def find_neighbor_file(sample_dir: Path, suffixes: Iterable[str]) -> Optional[Path]:
    for suf in suffixes:
        hits = sorted(sample_dir.glob(f"*_mtDNA_contig{suf}"))
        if hits:
            return hits[0]

    for suf in suffixes:
        hits = sorted(sample_dir.glob(f"*{suf}"))
        if hits:
            return hits[0]

    return None


def load_seq_from_any(seq_path: Path) -> tuple[str, Seq]:
    fmt = "genbank" if seq_path.suffix.lower() in {".gb", ".gbk", ".gbff"} else "fasta"
    rec = next(SeqIO.parse(str(seq_path), fmt))
    return rec.id, rec.seq


def parse_genbank_features(gb_path: Path) -> list[Feature]:
    rec = next(SeqIO.parse(str(gb_path), "genbank"))
    out: list[Feature] = []

    for feat in rec.features:
        if feat.type not in {"CDS", "tRNA", "rRNA"}:
            continue

        label = ""
        for key in ("gene", "product", "Name", "locus_tag", "label"):
            if key in feat.qualifiers and feat.qualifiers[key]:
                label = feat.qualifiers[key][0]
                break

        out.append(
            Feature(
                start=int(feat.location.start),
                end=int(feat.location.end),
                strand="+" if feat.location.strand != -1 else "-",
                ftype=feat.type,
                label=label,
            )
        )

    return out


def parse_gff_attributes(attr_text: str) -> dict[str, str]:
    out: dict[str, str] = {}

    for item in attr_text.strip().split(";"):
        if not item:
            continue

        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        elif " " in item:
            k, v = item.split(" ", 1)
            out[k] = v.strip('"')

    return out


def parse_gff_features(gff_path: Path) -> list[Feature]:
    out: list[Feature] = []

    with open(gff_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            _seqid, _source, ftype, start, end, _score, strand, _phase, attrs = parts

            if ftype not in {"CDS", "tRNA", "rRNA"}:
                continue

            ad = parse_gff_attributes(attrs)

            label = ""
            for key in ("gene", "product", "Name", "locus_tag", "ID"):
                if ad.get(key):
                    label = ad[key]
                    break

            out.append(
                Feature(
                    start=int(start) - 1,
                    end=int(end),
                    strand=strand if strand in {"+", "-"} else "+",
                    ftype=ftype,
                    label=label,
                )
            )

    return out


def load_genome_record(roots: list[Path], display_label: str) -> GenomeRecord:
    sample_code = extract_sample_code(display_label)
    sample_dir = get_best_sample_dir(roots, display_label)
    root = next((r for r in roots if str(sample_dir).startswith(str(r))), sample_dir.parents[2])

    gb_path = find_neighbor_file(sample_dir, [".gb", ".gbk", ".gbff"])
    gff_path = find_neighbor_file(sample_dir, [".gff", ".gff3"])
    fasta_path = find_neighbor_file(sample_dir, [".fasta", ".fa", ".fna"])

    if gb_path is None and fasta_path is None:
        raise FileNotFoundError(f"No GenBank or FASTA file found in {sample_dir}")

    if gb_path is not None:
        seq_id, seq = load_seq_from_any(gb_path)
        features = parse_genbank_features(gb_path)
        seq_path = gb_path
    else:
        if fasta_path is None:
            raise FileNotFoundError(f"No FASTA file found in {sample_dir}")

        seq_id, seq = load_seq_from_any(fasta_path)
        features = parse_gff_features(gff_path) if gff_path else []
        seq_path = fasta_path

    return GenomeRecord(
        display_label=display_label,
        sample_code=sample_code,
        root=root,
        sample_dir=sample_dir,
        seq_path=seq_path,
        seq_id=seq_id,
        seq=seq,
        gb_path=gb_path,
        gff_path=gff_path,
        features=features,
    )


# BLASTN alignment
def require_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found on PATH: {name}")


def write_fasta(path: Path, seq_id: str, seq: Seq) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f">{seq_id}\n")
        s = str(seq)
        for i in range(0, len(s), 80):
            fh.write(s[i:i + 80] + "\n")


def merge_intervals(intervals: list[tuple[int, int]], max_gap: int = 30) -> list[tuple[int, int]]:
    if not intervals:
        return []

    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged: list[list[int]] = [[intervals[0][0], intervals[0][1]]]

    for start, end in intervals[1:]:
        cur = merged[-1]
        if start <= cur[1] + max_gap:
            cur[1] = max(cur[1], end)
        else:
            merged.append([start, end])

    return [(s, e) for s, e in merged]


def run_blastn_nucleotide_alignment(
    reference: GenomeRecord,
    query: GenomeRecord,
    min_identity: float = 70.0,
    min_length: int = 80,
    merge_gap: int = 30,
) -> list[tuple[int, int]]:
    """Align one mitogenome to the reference and return reference intervals."""

    require_executable("makeblastdb")
    require_executable("blastn")

    with tempfile.TemporaryDirectory(prefix="mosq_blastn_nt_") as tmpdir:
        tmp = Path(tmpdir)

        ref_fa = tmp / "reference_nt.fasta"
        qry_fa = tmp / "query_nt.fasta"
        db_prefix = tmp / "reference_nt_db"
        out_tsv = tmp / "blastn_hits.tsv"

        write_fasta(ref_fa, reference.seq_id, reference.seq)
        write_fasta(qry_fa, query.seq_id, query.seq)

        subprocess.run(
            [
                "makeblastdb",
                "-in", str(ref_fa),
                "-dbtype", "nucl",
                "-out", str(db_prefix),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        subprocess.run(
            [
                "blastn",
                "-query", str(qry_fa),
                "-db", str(db_prefix),
                "-task", "blastn",
                "-dust", "no",
                "-evalue", "1e-20",
                "-max_target_seqs", "10000",
                "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send",
                "-out", str(out_tsv),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        intervals: list[tuple[int, int]] = []

        if out_tsv.exists() and out_tsv.stat().st_size > 0:
            with open(out_tsv, "r", encoding="utf-8") as fh:
                for line in fh:
                    (
                        _qseqid,
                        _sseqid,
                        pident,
                        length,
                        _qstart,
                        _qend,
                        sstart,
                        send,
                    ) = line.rstrip("\n").split("\t")

                    pident = float(pident)
                    length = int(length)
                    sstart = int(sstart)
                    send = int(send)

                    if pident < min_identity:
                        continue

                    if length < min_length:
                        continue

                    lo = min(sstart, send) - 1
                    hi = max(sstart, send)

                    intervals.append((lo, hi))

    return merge_intervals(intervals, max_gap=merge_gap)


# Circular plot
def pos_to_theta(pos: float, genome_len: int) -> float:
    return 90.0 - (pos / genome_len) * 360.0


def interval_to_thetas(start: int, end: int, genome_len: int) -> tuple[float, float]:
    return pos_to_theta(end, genome_len), pos_to_theta(start, genome_len)


def draw_ring_interval(
    ax,
    start: int,
    end: int,
    genome_len: int,
    r_outer: float,
    width: float,
    **kwargs,
):
    theta1, theta2 = interval_to_thetas(start, end, genome_len)
    ax.add_patch(Wedge((0, 0), r_outer, theta1, theta2, width=width, **kwargs))


def draw_full_ring(ax, genome_len: int, r_outer: float, width: float, **kwargs):
    draw_ring_interval(ax, 0, genome_len, genome_len, r_outer, width, **kwargs)


def midpoint(start: int, end: int) -> float:
    return start + (end - start) / 2.0


def polar_xy(radius: float, theta_deg: float) -> tuple[float, float]:
    theta_rad = math.radians(theta_deg)
    return radius * math.cos(theta_rad), radius * math.sin(theta_rad)


def text_on_circle(
    ax,
    text: str,
    pos: float,
    genome_len: int,
    radius: float,
    fontsize: float = 8.0,
    color: str = "black",
    fontstyle: Optional[str] = None,
):
    theta_deg = pos_to_theta(pos, genome_len)
    x, y = polar_xy(radius, theta_deg)

    rotation = theta_deg - 90.0
    ha = "left"

    if x < 0:
        rotation += 180.0
        ha = "right"

    ax.text(
        x,
        y,
        text,
        fontsize=fontsize,
        color=color,
        rotation=rotation,
        rotation_mode="anchor",
        ha=ha,
        va="center",
        fontstyle=fontstyle,
        clip_on=False,
    )


def clean_gene_label(label: str) -> str:
    label = label.strip()
    low = label.lower()

    if low.startswith("trn") and len(label) >= 4:
        aa = label[3:]
        aa = aa.replace("-", "")
        label = f"tRNA-{aa}"

    replacements = {
        "rrns": "s-rRNA",
        "rrnl": "l-rRNA",
        "12s ribosomal rna": "s-rRNA",
        "16s ribosomal rna": "l-rRNA",
    }

    if low in replacements:
        return replacements[low]

    label = label.replace("TRNA-", "tRNA-")
    return label


def estimate_label_sep_deg(text: str, fontsize: float) -> float:
    return max(5.2, min(18.0, 2.4 + len(text) * 0.42 + fontsize * 0.22))


def assign_label_lanes(labels: list[dict], max_lanes: int = 6) -> list[dict]:
    if not labels:
        return labels

    labels = sorted(labels, key=lambda d: d["angle360"])
    lane_last_angle = [-1e9] * max_lanes

    for lab in labels:
        min_sep = lab["min_sep_deg"]
        placed = False

        for lane in range(max_lanes):
            if lab["angle360"] - lane_last_angle[lane] >= min_sep:
                lab["lane"] = lane
                lane_last_angle[lane] = lab["angle360"]
                placed = True
                break

        if not placed:
            lab["lane"] = max_lanes - 1
            lane_last_angle[max_lanes - 1] = lab["angle360"]

    changed = True

    while changed:
        changed = False
        by_lane: dict[int, list[dict]] = {}

        for lab in labels:
            by_lane.setdefault(lab["lane"], []).append(lab)

        for lane, items in by_lane.items():
            items = sorted(items, key=lambda d: d["angle360"])

            if len(items) < 2:
                continue

            first = items[0]
            last = items[-1]
            wrap_gap = first["angle360"] + 360.0 - last["angle360"]
            need = max(first["min_sep_deg"], last["min_sep_deg"])

            if wrap_gap < need and lane < max_lanes - 1:
                last["lane"] += 1
                changed = True
                break

    return labels


def draw_label_connector(
    ax,
    pos: float,
    genome_len: int,
    r_from: float,
    r_to: float,
    lw: float = 0.55,
):
    theta_deg = pos_to_theta(pos, genome_len)
    x1, y1 = polar_xy(r_from, theta_deg)
    x2, y2 = polar_xy(r_to, theta_deg)

    ax.plot(
        [x1, x2],
        [y1, y2],
        color=LEADER_LINE_COLOR,
        lw=lw,
        alpha=0.85,
        solid_capstyle="round",
    )


def draw_reference_gene_ring(
    ax,
    reference: GenomeRecord,
    r_outer: float,
    width: float,
    label_radius: float,
    label_fontsize: float = 6.1,
    max_label_lanes: int = 6,
    lane_step: float = 0.034,
):
    features = sorted(reference.features or [], key=lambda f: (f.start, f.end))

    for feat in features:
        if feat.end <= feat.start:
            continue

        color = FEATURE_COLORS.get(feat.ftype, "#000000")

        draw_ring_interval(
            ax,
            feat.start,
            feat.end,
            reference.length,
            r_outer=r_outer,
            width=width,
            facecolor=color,
            edgecolor="white",
            linewidth=0.35,
        )

    label_items: list[dict] = []

    for feat in features:
        if feat.end <= feat.start or not feat.label:
            continue

        txt = clean_gene_label(feat.label)
        pos = midpoint(feat.start, feat.end)
        angle360 = (pos / reference.length) * 360.0

        fs = label_fontsize - 0.2 if feat.ftype in {"tRNA", "rRNA"} else label_fontsize

        label_items.append(
            {
                "text": txt,
                "pos": pos,
                "angle360": angle360,
                "fontsize": fs,
                "min_sep_deg": estimate_label_sep_deg(txt, fs),
                "lane": 0,
            }
        )

    label_items = assign_label_lanes(label_items, max_lanes=max_label_lanes)

    for lab in label_items:
        lane_radius = label_radius + lab["lane"] * lane_step

        draw_label_connector(
            ax,
            lab["pos"],
            reference.length,
            r_outer + 0.004,
            lane_radius - 0.015,
        )

        text_on_circle(
            ax,
            lab["text"],
            lab["pos"],
            reference.length,
            radius=lane_radius,
            fontsize=lab["fontsize"],
            color=GENE_LABEL_COLOR,
        )


def draw_inner_ticks(
    ax,
    genome_len: int,
    tick_every: int,
    radius: float,
    tick_len: float,
    fontsize: float = 7.0,
):
    for pos in range(0, genome_len + 1, tick_every):
        theta_deg = pos_to_theta(pos, genome_len)

        x1, y1 = polar_xy(radius, theta_deg)
        x2, y2 = polar_xy(radius - tick_len, theta_deg)

        ax.plot([x1, x2], [y1, y2], color=INNER_TICK_COLOR, lw=0.5)

        if pos == 0:
            continue

        text_on_circle(
            ax,
            f"{pos // 1000} kbp",
            pos,
            genome_len,
            radius - tick_len - 0.045,
            fontsize=fontsize,
            color=INNER_TICK_COLOR,
        )


def wrap_center_label(label: str, width: int = 20) -> str:
    species = re.sub(r"\s*\(F\d+_\d+\)$", "", label)
    words = species.split()

    lines = []
    cur = []
    cur_len = 0

    for w in words:
        if cur and cur_len + 1 + len(w) > width:
            lines.append(" ".join(cur))
            cur = [w]
            cur_len = len(w)
        else:
            cur.append(w)
            cur_len += len(w) + (1 if len(cur) > 1 else 0)

    if cur:
        lines.append(" ".join(cur))

    return "\n".join(lines)


def make_plot(
    reference: GenomeRecord,
    genomes: list[GenomeRecord],
    alignments: dict[str, list[tuple[int, int]]],
    output_png: Path,
    output_svg: Path,
    figure_size: float = 11.6,
    dpi: int = 600,
    tick_every: int = 1000,
    gene_track_width: float = 0.040,
    ring_width: float = 0.026,
    ring_gap: float = 0.0042,
):
    fig, ax = plt.subplots(figsize=(figure_size, figure_size))
    ax.set_aspect("equal")
    ax.axis("off")

    outer_radius = 1.0
    gene_outer = outer_radius
    gene_label_radius = gene_outer + 0.065

    draw_reference_gene_ring(
        ax,
        reference,
        r_outer=gene_outer,
        width=gene_track_width,
        label_radius=gene_label_radius,
        label_fontsize=5.8,
        max_label_lanes=10,
        lane_step=0.048,
    )

    start_outer = gene_outer - gene_track_width - 0.012

    for i, genome in enumerate(genomes):
        r_outer = start_outer - i * (ring_width + ring_gap)

        draw_full_ring(
            ax,
            genome_len=reference.length,
            r_outer=r_outer,
            width=ring_width,
            facecolor=BACKGROUND_TRACK_COLOR,
            edgecolor="white",
            linewidth=0.15,
        )

        for start, end in alignments.get(genome.display_label, []):
            draw_ring_interval(
                ax,
                start=start,
                end=end,
                genome_len=reference.length,
                r_outer=r_outer,
                width=ring_width,
                facecolor=SELECTED_SAMPLE_COLORS[genome.display_label],
                edgecolor="none",
            )

    inner_radius = (
        start_outer
        - (len(genomes) - 1) * (ring_width + ring_gap)
        - ring_width
        - 0.018
    )

    draw_inner_ticks(
        ax,
        reference.length,
        tick_every=tick_every,
        radius=inner_radius + 0.012,
        tick_len=0.022,
        fontsize=7.0,
    )

    center_species = wrap_center_label(reference.display_label, width=20)
    center_species = re.sub(r"\s*\(F\d+_\d+\)$", "", center_species)

    ax.text(
        0,
        0,
        f"{center_species}\n{reference.length:,} bp",
        ha="center",
        va="center",
        fontsize=10.5,
        color=CENTER_TEXT_COLOR,
        fontstyle="italic",
    )

    handles = []
    labels = []

    # Genome legend
    for idx, genome in enumerate(genomes, start=1):
        handles.append(
            Rectangle(
                (0, 0),
                1,
                1,
                facecolor=SELECTED_SAMPLE_COLORS[genome.display_label],
                edgecolor="none",
            )
        )
        labels.append(f"{idx}. {genome.display_label}")

    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=4,
        frameon=False,
        fontsize=7.1,
        handlelength=0.95,
        handletextpad=0.34,
        columnspacing=1.1,
        labelspacing=0.55,
    )

    # Gene legend
    feature_handles = [
        Rectangle((0, 0), 1, 1, facecolor=FEATURE_COLORS["CDS"], edgecolor="none"),
        Rectangle((0, 0), 1, 1, facecolor=FEATURE_COLORS["tRNA"], edgecolor="none"),
        Rectangle((0, 0), 1, 1, facecolor=FEATURE_COLORS["rRNA"], edgecolor="none"),
    ]

    feature_labels = [
        "PCGs",
        "tRNAs",
        "rRNAs",
    ]

    ax.legend(
        feature_handles,
        feature_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.04),
        ncol=3,
        frameon=False,
        fontsize=8.0,
        handlelength=1.2,
        handletextpad=0.4,
        columnspacing=1.0,
    )

    lim = gene_label_radius + 10 * 0.048 + 0.13
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    fig.subplots_adjust(top=0.955, bottom=0.19, left=0.03, right=0.97)

    output_png.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_png, dpi=dpi, bbox_inches="tight")
    fig.savefig(output_svg, bbox_inches="tight")

    plt.close(fig)


# CLI
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Comparative BRIG-like plot for selected mosquito mitogenomes using nucleotide BLASTN alignment"
    )

    p.add_argument(
        "--roots",
        nargs="+",
        default=[
            "/path/to/FINAL_mosquitos",
        ],
        help="Project roots; each root must contain a data directory",
    )

    p.add_argument(
        "--outdir",
        default="./comparative_brig_plot",
        help="Output directory",
    )

    p.add_argument(
        "--prefix",
        default="selected_mosquitoes_comparative_nt_alignment",
        help="Output filename prefix",
    )

    p.add_argument(
        "--min-identity",
        type=float,
        default=70.0,
        help="Minimum BLASTN identity percentage",
    )

    p.add_argument(
        "--min-length",
        type=int,
        default=80,
        help="Minimum BLASTN alignment length",
    )

    p.add_argument(
        "--merge-gap",
        type=int,
        default=30,
        help="Merge blocks separated by at most this many nucleotides",
    )

    p.add_argument(
        "--tick-every",
        type=int,
        default=1000,
        help="Spacing between internal tick marks in base pairs",
    )

    return p.parse_args()


def main() -> None:
    args = parse_args()

    roots = [Path(x).expanduser().resolve() for x in args.roots]
    outdir = Path(args.outdir).expanduser().resolve()

    out_png = outdir / f"{args.prefix}.png"
    out_svg = outdir / f"{args.prefix}.svg"

    print("[INFO] Loading reference...")
    reference = load_genome_record(roots, REFERENCE_LABEL)

    if not reference.features:
        raise RuntimeError(
            f"The reference {REFERENCE_LABEL} has no readable CDS, tRNA, or rRNA features. "
            f"Provide a valid GenBank or GFF file in {reference.sample_dir}."
        )

    print(f"[INFO] Reference: {REFERENCE_LABEL}")
    print(f"[INFO] Reference length: {reference.length:,} bp")
    print(f"[INFO] Reference file: {reference.seq_path}")

    print("[INFO] Loading selected genomes...")
    genomes: list[GenomeRecord] = []

    for label in COMPARISON_LABELS:
        rec = load_genome_record(roots, label)
        genomes.append(rec)
        print(f"  - {label} -> {rec.sample_dir}")

    print("[INFO] Running BLASTN nucleotide alignments against the reference...")

    alignments: dict[str, list[tuple[int, int]]] = {}

    for genome in genomes:
        intervals = run_blastn_nucleotide_alignment(
            reference=reference,
            query=genome,
            min_identity=args.min_identity,
            min_length=args.min_length,
            merge_gap=args.merge_gap,
        )

        alignments[genome.display_label] = intervals

        cov = sum(e - s for s, e in intervals)
        pct = 100.0 * cov / reference.length if reference.length else 0.0

        print(
            f"  - {genome.display_label}: "
            f"{len(intervals)} nucleotide blocks, reference coverage about {pct:.2f}%"
        )

    print("[INFO] Drawing figure...")

    make_plot(
        reference=reference,
        genomes=genomes,
        alignments=alignments,
        output_png=out_png,
        output_svg=out_svg,
        tick_every=args.tick_every,
    )

    print(f"[OK] PNG: {out_png}")
    print(f"[OK] SVG: {out_svg}")


if __name__ == "__main__":
    main()
