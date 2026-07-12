#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Generate Circos plots from corrected MitoFinder GFF and FASTA outputs.

from pycirclize import Circos, config
from pycirclize.parser import Gff
from Bio import SeqIO
import numpy as np
from pathlib import Path
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import re
import textwrap
import sys
import argparse


# Plot configuration.
feature_cfg = {
    "CDS":  {"rlim": (95, 100), "color": "#D55E00", "style": "arrow"},
    "tRNA": {"rlim": (90, 95),  "color": "#009E73", "style": "box"},
    "rRNA": {"rlim": (90, 95),  "color": "#0072B2", "style": "box"},
}

HEAVY_IS_PLUS = True

CDS_LABEL_R = (101.0, 109.0)
RNA_LABEL_R = (96.0, 104.0)

try:
    config.ann_adjust.enable = True
except Exception:
    pass

GC_TRACK_RLIM = (75, 85)
GC_POS_COLOR  = "olive"
GC_NEG_COLOR  = "purple"
GC_WINDOW     = 15
GC_STEP       = 20
GC_ALPHA      = 2
GC_SKEW_FLIP  = False
IS_CIRCULAR   = True
GC_VMIN, GC_VMAX = -0.7, 0.7
GC_EPS        = 1e-9

EDGE_COLOR   = "none"
AXIS_FILL    = "#F7F7F7"
TICK_EC      = "grey"
TICK_LW      = 0.6
LABEL_SIZE   = 7

INCLUDE_HYPOTHETICAL = False


# Helpers.
def get_label(qual):
    for key in ("gene", "product", "Name", "locus_tag", "ID"):
        if key in qual and qual[key]:
            return qual[key][0]
    return ""


def short_label(text, maxlen=22):
    return text if len(text) <= maxlen else text[:maxlen] + "..."


def get_feat_strand(feat):
    loc = getattr(feat, "location", None)
    return getattr(loc, "strand", None)


def is_heavy_strand(feat):
    s = get_feat_strand(feat)
    if s not in (1, -1):
        return False
    return (s == 1) if HEAVY_IS_PLUS else (s == -1)


def compute_gc_skew_circular(seq, win=GC_WINDOW, step=GC_STEP, alpha=GC_ALPHA, flip=GC_SKEW_FLIP):
    s = str(seq).upper()
    n = len(s)
    if n == 0:
        return np.array([0]), np.array([0.0])

    win = max(1, min(win, max(1, n - 1)))
    step = max(1, step)
    s2 = s + s[:win - 1]

    pos, skew = [], []

    for i in range(0, n, step):
        w = s2[i:i + win]
        g, c = w.count("G"), w.count("C")
        g_s, c_s = g + alpha, c + alpha
        val = (g_s - c_s) / (g_s + c_s)

        if flip:
            val = -val

        pos.append((i + win // 2) % n)
        skew.append(val)

    idx = np.argsort(pos)
    return np.array(pos)[idx], np.array(skew)[idx]


def compute_gc_skew_linear(seq, win=GC_WINDOW, step=GC_STEP, alpha=GC_ALPHA, flip=GC_SKEW_FLIP):
    s = str(seq).upper()
    n = len(s)
    if n == 0:
        return np.array([0]), np.array([0.0])

    win = max(1, min(win, n))
    step = max(1, step)

    pos, skew = [], []

    for i in range(0, max(1, n - win + 1), step):
        w = s[i:i + win]
        g, c = w.count("G"), w.count("C")
        g_s, c_s = g + alpha, c + alpha
        val = (g_s - c_s) / (g_s + c_s)

        if flip:
            val = -val

        pos.append(min(i + win // 2, n - 1))
        skew.append(val)

    if not pos:
        pos = [n // 2]
        skew = [0.0]

    return np.array(pos), np.array(skew)


def load_fasta_as_dict(fasta_path: Path):
    fasta_path = Path(fasta_path)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file was not found: {fasta_path}")

    return {rec.id: rec.seq for rec in SeqIO.parse(str(fasta_path), "fasta")}


def genus_from_path(gff_path: Path, data_dir: Path) -> str:
    """
    Expected layout:
      results/annotation/corrections/medaka_mitofinder_rrnS/culex/F1_83-Culex_quinquefasciatus/*.gff

    Returns:
      culex
    """
    try:
        rel = gff_path.resolve().relative_to(data_dir.resolve())
        return rel.parts[0]
    except Exception:
        return "unknown"


def sample_name_from_path(gff_path: Path, data_dir: Path) -> str:
    """
    Expected layout:
      results/annotation/corrections/medaka_mitofinder_rrnS/culex/F1_83-Culex_quinquefasciatus/*.gff

    Returns:
      F1_83-Culex_quinquefasciatus
    """
    try:
        rel = gff_path.resolve().relative_to(data_dir.resolve())

        if len(rel.parts) >= 2:
            return rel.parts[1]

    except Exception:
        pass

    return gff_path.parent.name


def sample_id_from_sample_dir(sample_dir_name: str) -> str:
    """
    F1_83-Culex_quinquefasciatus -> F1_83
    F3_95-Culex_iolambdis       -> F3_95
    F4_93-Culex_sp.             -> F4_93
    """
    if "-" in sample_dir_name:
        return sample_dir_name.split("-", 1)[0]

    m = re.search(r"(F\d+[_-]?\d+)", sample_dir_name, flags=re.IGNORECASE)
    if m:
        return m.group(1).replace("-", "_").upper()

    return sample_dir_name


def mosquito_name_from_sample(sample_dir_name: str) -> str:
    """
    F1_83-Culex_quinquefasciatus -> Culex quinquefasciatus
    F4_93-Culex_sp.              -> Culex sp.
    """
    if "-" in sample_dir_name:
        mosq = sample_dir_name.split("-", 1)[1]
    else:
        mosq = sample_dir_name

    mosq = mosq.replace("_", " ")
    mosq = re.sub(r"\s+", " ", mosq).strip()

    return mosq


def barcode_from_filename(filename: str) -> str:
    """
    Extract the barcode token from the result filename.
    F1_barcode83_rrnS_mtDNA_contig.gff -> barcode83
    barcode83_rrnS_mtDNA_contig.gff    -> barcode83
    """
    name = Path(filename).name

    m = re.search(r"(barcode\d+)", name, flags=re.IGNORECASE)

    if m:
        return m.group(1).lower()

    return Path(filename).stem


def sample_id_from_filename(filename: str) -> str:
    """
    F1_barcode83_rrnS_mtDNA_contig.gff -> F1_83
    F3_barcode95_rrnS_mtDNA_contig.gff -> F3_95
    """
    name = Path(filename).name

    m = re.search(
        r"(?P<run>F\d+)[_-](?P<barcode>barcode(?P<num>\d+))",
        name,
        flags=re.IGNORECASE
    )

    if m:
        run = m.group("run").upper()
        num = m.group("num")
        return f"{run}_{num}"

    m = re.search(r"(barcode)(\d+)", name, flags=re.IGNORECASE)

    if m:
        return f"barcode{m.group(2)}"

    return Path(filename).stem


def italics_mathtext(text: str) -> str:
    t = text.replace("\\", r"\\").replace("_", r"\_")
    t = t.replace(" ", r"\;")
    return rf"$\it{{{t}}}$"


def format_center_title(name: str, wrap_width: int = 22):
    lines = textwrap.wrap(name, width=wrap_width, break_long_words=False)

    if not lines:
        lines = [name]

    n_lines = len(lines)
    max_line_len = max(len(x) for x in lines)

    fs = 15

    if n_lines == 2:
        fs = 14
    elif n_lines == 3:
        fs = 13
    elif n_lines >= 4:
        fs = 12

    if max_line_len >= 26:
        fs -= 2
    elif max_line_len >= 22:
        fs -= 1

    fs = max(9, fs)

    return "\n".join(lines), fs


def total_nucleotides(seq_dict, seqid2size) -> int:
    total = 0

    for seqid in seqid2size.keys():
        if seqid in seq_dict and seq_dict[seqid] is not None:
            total += len(seq_dict[seqid])
        else:
            total += int(seqid2size[seqid])

    return total


# =====================================================
# Plot
# =====================================================
def plot_one_barcode(
    gff_path: Path,
    fasta_path: Path,
    out_prefix: Path,
    sample_dir_name: str,
    wrap_width: int = 22
):
    sample_id = sample_id_from_sample_dir(sample_dir_name)
    barcode = barcode_from_filename(gff_path.name)
    mosq_name = mosquito_name_from_sample(sample_dir_name)

    print(f"[INFO] Processing: {sample_id} | {barcode} | {mosq_name}")

    gff = Gff(str(gff_path))
    seqid2size = gff.get_seqid2size()
    seq_dict = load_fasta_as_dict(fasta_path)

    if len(seqid2size) == 1 and len(seq_dict) == 1:
        gff_id = list(seqid2size.keys())[0]
        fa_id = list(seq_dict.keys())[0]

        if gff_id != fa_id:
            seq_dict = {gff_id: seq_dict[fa_id]}

    seqid2features_all = {
        ftype: gff.get_seqid2features(feature_type=ftype)
        for ftype in feature_cfg
    }

    title_text, title_size = format_center_title(mosq_name, wrap_width=wrap_width)
    title_text = "\n".join(
        italics_mathtext(line)
        for line in title_text.splitlines()
    )

    n_nt = total_nucleotides(seq_dict, seqid2size)
    title_text = f"{title_text}\n(length = {n_nt} bp)"
    title_size = max(8, title_size - 1)

    circos = Circos(sectors=seqid2size, space=0 if len(seqid2size) == 1 else 2)
    circos.text(title_text, size=title_size)

    for sector in circos.sectors:
        seq = seq_dict.get(sector.name)

        if seq is not None:
            gc_skew_track = sector.add_track(GC_TRACK_RLIM)
            gc_skew_track.axis(fc=AXIS_FILL, ec="none")

            if IS_CIRCULAR:
                label_pos_list, gc_skews = compute_gc_skew_circular(
                    seq,
                    win=GC_WINDOW,
                    step=GC_STEP,
                    alpha=GC_ALPHA,
                    flip=GC_SKEW_FLIP
                )
            else:
                label_pos_list, gc_skews = compute_gc_skew_linear(
                    seq,
                    win=GC_WINDOW,
                    step=GC_STEP,
                    alpha=GC_ALPHA,
                    flip=GC_SKEW_FLIP
                )

            gc_skews = np.clip(gc_skews, GC_VMIN + GC_EPS, GC_VMAX - GC_EPS)

            positive_gc_skews = np.where(gc_skews > 0, gc_skews, 0)
            negative_gc_skews = np.where(gc_skews < 0, gc_skews, 0)

            gc_skew_track.fill_between(
                label_pos_list,
                positive_gc_skews,
                0,
                vmin=GC_VMIN,
                vmax=GC_VMAX,
                color=GC_POS_COLOR
            )

            gc_skew_track.fill_between(
                label_pos_list,
                negative_gc_skews,
                0,
                vmin=GC_VMIN,
                vmax=GC_VMAX,
                color=GC_NEG_COLOR
            )

            gc_skew_track.xticks_by_interval(
                interval=1000,
                outer=False,
                label_formatter=lambda v: f"{v / 1000:.1f} Kb",
                label_orientation="vertical",
                line_kws=dict(ec=TICK_EC, lw=TICK_LW),
            )

        heavy_track = sector.add_track((95, 100))
        heavy_track.axis(fc=AXIS_FILL, ec="none")

        light_track = sector.add_track((90, 95))
        light_track.axis(fc=AXIS_FILL, ec="none")

        for ftype in ("CDS", "tRNA", "rRNA"):
            feats = seqid2features_all[ftype].get(sector.name, [])
            cfg = feature_cfg[ftype]

            for feat in feats:
                label = get_label(feat.qualifiers)

                if (
                    not INCLUDE_HYPOTHETICAL
                    and label
                    and label.lower().startswith("hypothetical")
                ):
                    continue

                strand = get_feat_strand(feat)

                if strand not in (1, -1):
                    continue

                target_track = heavy_track if is_heavy_strand(feat) else light_track

                target_track.genomic_features(
                    feat,
                    plotstyle=cfg["style"],
                    fc=cfg["color"],
                    ec=EDGE_COLOR
                )

                if label:
                    start, end = int(feat.location.start), int(feat.location.end)
                    pos = (start + end) // 2

                    label_r = CDS_LABEL_R if ftype == "CDS" else RNA_LABEL_R

                    target_track.annotate(
                        pos,
                        short_label(label, maxlen=22),
                        min_r=label_r[0],
                        max_r=label_r[1],
                        label_size=LABEL_SIZE,
                        line_kws=dict(lw=0.6, alpha=0.85),
                        text_kws=dict(va="center"),
                    )

    handles = [
        Patch(facecolor=feature_cfg["CDS"]["color"],  label="CDS"),
        Patch(facecolor=feature_cfg["tRNA"]["color"], label="tRNA"),
        Patch(facecolor=feature_cfg["rRNA"]["color"], label="rRNA"),
        Patch(facecolor=GC_POS_COLOR, label="GC skew > 0"),
        Patch(facecolor=GC_NEG_COLOR, label="GC skew < 0"),
    ]

    fig = circos.plotfig(figsize=(9, 9))

    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02)
    )

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(str(out_prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(out_prefix) + ".png", dpi=600, bbox_inches="tight")

    plt.close(fig)

    print(f"[OK] Saved: {out_prefix}.svg")
    print(f"[OK] Saved: {out_prefix}.png")


# Find GFF/FASTA pairs in the corrected result layout.
def find_gff_fasta_pairs(data_dir: Path, gff_pattern="*_mtDNA_contig.gff"):
    """
    Search recursively in:

      results/annotation/corrections/medaka_mitofinder_rrnS/
        culex/
          F1_83-Culex_quinquefasciatus/
            F1_barcode83_rrnS_mtDNA_contig.gff
            F1_barcode83_rrnS_mtDNA_contig.fasta

    Returns:
      genus, sample_dir, sample_id, barcode, gff_path, fasta_path
    """
    pairs = []

    for gff_path in sorted(data_dir.rglob(gff_pattern)):
        gff_path = gff_path.resolve()
        fasta_path = gff_path.with_suffix(".fasta")

        if not fasta_path.exists():
            print(f"[WARN] FASTA file not found for: {gff_path}", file=sys.stderr)
            continue

        genus = genus_from_path(gff_path, data_dir)
        sample_dir = sample_name_from_path(gff_path, data_dir)

        sample_id_folder = sample_id_from_sample_dir(sample_dir)
        sample_id_file = sample_id_from_filename(gff_path.name)

        barcode = barcode_from_filename(gff_path.name)

        if sample_id_folder != sample_id_file:
            print(
                f"[WARN] Folder and file sample IDs do not match: "
                f"folder={sample_id_folder} | file={sample_id_file} | {gff_path}",
                file=sys.stderr
            )

        sample_id = sample_id_folder

        pairs.append(
            (
                genus,
                sample_dir,
                sample_id,
                barcode,
                gff_path,
                fasta_path
            )
        )

    return pairs


# =====================================================
# Main
# =====================================================
def main():
    script_dir = Path(__file__).resolve().parent
    workspace_root = script_dir.parent

    p = argparse.ArgumentParser(
        description=(
            "Generate Circos plots for mitochondrial genomes from corrected MitoFinder outputs."
        )
    )

    p.add_argument(
        "--root",
        default=str(workspace_root),
        help="Workspace root. Default: parent directory of this script."
    )

    p.add_argument(
        "--data-subdir",
        default="results/annotation/corrections/medaka_mitofinder_rrnS",
        help="Relative path to the corrected MitoFinder result tree."
    )

    p.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Default: <root>/results/annotation/plots/circos"
    )

    p.add_argument(
        "--gff-pattern",
        default="*_mtDNA_contig.gff",
        help="GFF filename pattern. Default: *_mtDNA_contig.gff"
    )

    p.add_argument(
        "--wrap-width",
        type=int,
        default=22,
        help="Wrapping width for the center taxon label. Default: 22"
    )

    args = p.parse_args()

    root = Path(args.root).expanduser().resolve()
    data_dir = (root / args.data_subdir).resolve()

    if not data_dir.is_dir():
        raise SystemExit(f"[ERROR] Data directory does not exist: {data_dir}")

    if args.outdir:
        output_dir = Path(args.outdir).expanduser().resolve()
    else:
        output_dir = (root / "results" / "annotation" / "plots" / "circos").resolve()

    out_base = output_dir / "by_genus"

    pairs = find_gff_fasta_pairs(
        data_dir,
        gff_pattern=args.gff_pattern
    )

    if not pairs:
        raise SystemExit(f"[ERROR] No GFF/FASTA pairs were found in: {data_dir}")

    print(f"[INFO] Found {len(pairs)} mitogenomes.")

    for genus, sample_dir, sample_id, barcode, gff_path, fasta_path in pairs:
        out_prefix = out_base / genus / sample_dir / sample_id

        plot_one_barcode(
            gff_path=gff_path,
            fasta_path=fasta_path,
            out_prefix=out_prefix,
            sample_dir_name=sample_dir,
            wrap_width=args.wrap_width
        )


if __name__ == "__main__":
    main()
