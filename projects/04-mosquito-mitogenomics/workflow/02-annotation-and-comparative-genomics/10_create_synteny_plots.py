#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Generate per-genus and global synteny plots from corrected GenBank outputs.

import argparse
import sys
import re
import shutil
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib
from matplotlib.colors import Normalize

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank

from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqFeature import FeatureLocation, CompoundLocation


CDS_COLOR = "#D55E00"
TRNA_COLOR = "#009E73"
RRNA_COLOR = "#0072B2"


@dataclass(frozen=True)
class GeneHit:
    """A gene feature projected on a track/segment with sequence for identity calc."""
    gene_id: str
    ftype: str
    seqid: str
    start_1based: int
    end_1based: int
    strand: int
    seq: str


def parse_args():
    script_dir = Path(__file__).resolve().parent
    workspace_root = script_dir.parent

    p = argparse.ArgumentParser(
        description=(
            "Generate per-genus and global mitochondrial synteny plots from corrected GenBank outputs."
        )
    )

    p.add_argument(
        "--root",
        default=str(workspace_root),
        help="Workspace root. Default: parent directory of this script.",
    )
    p.add_argument(
        "--data-subdir",
        default="results/annotation/corrections/medaka_mitofinder_rrnS",
        help="Relative path to the corrected MitoFinder GenBank result tree.",
    )
    p.add_argument(
        "--pattern",
        default="*_mtDNA_contig*.gb",
        help="GenBank filename pattern (default: '*_mtDNA_contig*.gb')",
    )

    p.add_argument(
        "--ref-dir",
        default=str(workspace_root / "results" / "annotation" / "ref_gb"),
        help="Directory that contains reference mitogenomes.",
    )
    p.add_argument(
        "--ref-pattern",
        default="*.gb",
        help="GenBank pattern inside --ref-dir (default: '*.gb')",
    )

    p.add_argument(
        "--no-rotate-ref-to-rrnS",
        action="store_true",
        help="Disable reference rotation to rrnS end + 10 nt.",
    )
    p.add_argument(
        "--rrns-regex",
        default=r"(?:\brrnS\b|\b12S\b|12S[-\s]*rRNA|small[-\s]*subunit[-\s]*(?:ribosomal\s*)?RNA|s[-\s]*rRNA)",
        help="Regex used to detect rrnS/12S rRNA in feature qualifiers.",
    )

    p.add_argument("--identity-thr", type=float, default=70.0, help="Minimum identity percent (default 70).")
    p.add_argument("--length-thr-nt", type=int, default=70, help="Minimum nucleotide length (default 70).")

    p.add_argument("--barsize", type=int, default=2000, help="Scale bar length in bp (default 2000).")
    p.add_argument("--dpi", type=int, default=300, help="DPI PNG (default 300).")
    p.add_argument("--alpha", type=float, default=0.85, help="Alpha ribbons (default 0.85).")
    p.add_argument("--show-trna-rrna", action="store_true", help="Draw tRNA/rRNA features and gene links.")

    p.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Default: <root>/results/annotation/plots/synteny",
    )
    p.add_argument(
        "--prefix-all",
        default="all_synteny_nt",
        help="Filename prefix for the global output (default: all_synteny_nt)",
    )
    p.add_argument(
        "--prefix-suffix",
        default="synteny_nt",
        help="Per-genus filename suffix: <genus>_<suffix> (default: synteny_nt)",
    )
    return p.parse_args()


def list_genera(data_dir: Path) -> List[Path]:
    if not data_dir.is_dir():
        sys.exit(f"[ERROR] Data directory does not exist: {data_dir}")
    genera = sorted([p for p in data_dir.iterdir() if p.is_dir()])
    if not genera:
        sys.exit(f"[ERROR] No genus subdirectories were found in: {data_dir}")
    return genera


def find_gb_in_dir(base: Path, pattern: str) -> List[Path]:
    return sorted(base.rglob(pattern))


def sample_dir_from_gb(gb_path: Path, data_dir: Path, ref_dir: Optional[Path] = None) -> str:
    """
    For the corrected data tree:
      <data>/<genus>/<sample_dir>/<file.gb> -> sample_dir
    For the reference tree:
      <ref_dir>/<ACC>-<Taxon>/<ACC>.gb -> "<ACC>-<Taxon>"
    """
    gb_path = gb_path.resolve()
    data_dir = data_dir.resolve()
    ref_dir = ref_dir.resolve() if ref_dir else None

    if ref_dir and (ref_dir in gb_path.parents):
        try:
            rel = gb_path.relative_to(ref_dir)
            parts = rel.parts
            if len(parts) >= 1:
                return parts[0]
        except Exception:
            pass
        return gb_path.stem

    try:
        rel = gb_path.relative_to(data_dir)
        parts = rel.parts
        if len(parts) >= 2:
            return parts[1]
    except Exception:
        pass

    return gb_path.stem


ID_TAXON_RE = re.compile(r"^(?P<sid>F\d+[_-]\d+)[-_]+(?P<taxon>.+)$")
REF_ACC_TAXON_RE = re.compile(r"^(?P<acc>[A-Za-z]{1,4}\d+)[-_]+(?P<taxon>.+)$")


def _format_taxon_name(raw: str) -> str:
    if not raw:
        return raw

    taxon = raw.replace("_", " ").replace("-", " ").strip()
    taxon = re.sub(r"\s+", " ", taxon)

    taxon = re.sub(r"\bsp\b\.?", "sp.", taxon, flags=re.IGNORECASE)
    taxon = re.sub(r"\bspp\b\.?", "spp.", taxon, flags=re.IGNORECASE)
    taxon = re.sub(r"\bcf\b\.?", "cf.", taxon, flags=re.IGNORECASE)
    taxon = re.sub(r"\baff\b\.?", "aff.", taxon, flags=re.IGNORECASE)
    taxon = re.sub(r"\bnr\b\.?", "nr.", taxon, flags=re.IGNORECASE)

    parts = taxon.split(" ")
    if not parts:
        return taxon

    parts[0] = parts[0][0].upper() + parts[0][1:] if len(parts[0]) > 1 else parts[0].upper()

    for i in range(1, len(parts)):
        token = parts[i]
        if token.lower() in {"sp.", "spp.", "cf.", "aff.", "nr."}:
            parts[i] = token.lower()
        elif token.isalpha():
            parts[i] = token.lower()
        else:
            parts[i] = token

    return " ".join(parts)


def _italic_mathtext(text: str) -> str:
    if not text:
        return text
    s = text.replace("\\", r"\\")
    s = s.replace(" ", r"\ ")
    s = s.replace("_", r"\_")
    return rf"$\it{{{s}}}$"


def track_label_plain_from_gb(gb_path: Path, data_dir: Path, ref_dir: Optional[Path] = None) -> str:
    sample_dir = sample_dir_from_gb(gb_path, data_dir, ref_dir=ref_dir)
    if not sample_dir:
        return gb_path.stem

    m = ID_TAXON_RE.match(sample_dir)
    if m:
        sid = m.group("sid").replace("-", "_")
        taxon = _format_taxon_name(m.group("taxon"))
        return f"{taxon} ({sid})"

    m2 = REF_ACC_TAXON_RE.match(sample_dir)
    if m2:
        acc = m2.group("acc")
        taxon = _format_taxon_name(m2.group("taxon"))
        return f"{taxon} ({acc})"

    return sample_dir


def sort_key_synteny_custom(gb_path: Path, data_dir: Path, ref_dir: Optional[Path] = None):
    """
    Keep a stable Culex ordering while preserving IDs in the final labels.
      1) Culex iolambdis
      2) Culex sp.
      3) Culex quinquefasciatus
      4) other Culex

    Example labels remain:
      Culex sp. (F4_90)
      Culex quinquefasciatus (F1_83)
    """
    label = track_label_plain_from_gb(gb_path, data_dir, ref_dir=ref_dir)

    # Keep the ID in the display label, but sort by the taxon text only.
    taxon = re.sub(r"\s*\([^)]*\)\s*$", "", label).strip()
    taxon_low = taxon.lower()

    parts = taxon_low.split()
    genus = parts[0] if parts else taxon_low

    if genus == "culex":
        if taxon_low.startswith("culex iolambdis"):
            species_rank = 0
        elif taxon_low.startswith("culex sp"):
            species_rank = 1
        elif taxon_low.startswith("culex quinquefasciatus"):
            species_rank = 2
        else:
            species_rank = 3

        return (genus, species_rank, taxon_low, label.lower())

    return (genus, 0, taxon_low, label.lower())


def track_label_from_gb(gb_path: Path, data_dir: Path, ref_dir: Optional[Path] = None) -> str:
    sample_dir = sample_dir_from_gb(gb_path, data_dir, ref_dir=ref_dir)
    if not sample_dir:
        return gb_path.stem

    m = ID_TAXON_RE.match(sample_dir)
    if m:
        sid = m.group("sid").replace("-", "_")
        taxon_plain = _format_taxon_name(m.group("taxon"))
        taxon_ital = _italic_mathtext(taxon_plain)
        return f"{taxon_ital} ({sid})"

    m2 = REF_ACC_TAXON_RE.match(sample_dir)
    if m2:
        acc = m2.group("acc")
        taxon_plain = _format_taxon_name(m2.group("taxon"))
        taxon_ital = _italic_mathtext(taxon_plain)
        return f"{taxon_ital} ({acc})"

    return sample_dir


def _qual_text(qdict) -> str:
    if not qdict:
        return ""
    parts = []
    for v in qdict.values():
        if isinstance(v, list):
            parts.extend(map(str, v))
        else:
            parts.append(str(v))
    return " ".join(parts)


def find_rrnS_cut(rec, rrns_re: re.Pattern) -> Optional[int]:
    """
    Return 0-based cut position corresponding to rrnS end + 10 nt.
    Searches rRNA features first, then gene features.
    """
    L = len(rec.seq)
    if L <= 0:
        return None

    for ft in rec.features:
        if ft.type != "rRNA":
            continue
        txt = _qual_text(getattr(ft, "qualifiers", {}) or {})
        if rrns_re.search(txt):
            return (int(ft.location.end) + 10) % L

    for ft in rec.features:
        if ft.type != "gene":
            continue
        txt = _qual_text(getattr(ft, "qualifiers", {}) or {})
        if rrns_re.search(txt):
            return (int(ft.location.end) + 10) % L

    return None


def shift_location(loc, cut: int, L: int):
    strand = loc.strand
    if isinstance(loc, CompoundLocation):
        new_parts = []
        for part in loc.parts:
            shifted = shift_location(part, cut, L)
            if isinstance(shifted, CompoundLocation):
                new_parts.extend(list(shifted.parts))
            else:
                new_parts.append(shifted)
        return CompoundLocation(new_parts, operator="join") if len(new_parts) > 1 else new_parts[0]

    s = (int(loc.start) - cut) % L
    e = (int(loc.end) - cut) % L

    if s < e:
        return FeatureLocation(s, e, strand=strand)
    elif s > e:
        p1 = FeatureLocation(s, L, strand=strand)
        p2 = FeatureLocation(0, e, strand=strand)
        return CompoundLocation([p1, p2], operator="join")
    else:
        return FeatureLocation(0, 0, strand=strand)


def rotate_record_to_cut(rec, cut: int):
    """
    Rotate circular sequence so that new position 0 corresponds to old 'cut'.
    Adjust all features accordingly.
    """
    from copy import deepcopy

    L = len(rec.seq)
    if cut <= 0 or cut >= L:
        return rec

    new_rec = deepcopy(rec)
    new_rec.seq = rec.seq[cut:] + rec.seq[:cut]
    new_rec.features = []
    for ft in rec.features:
        nft = deepcopy(ft)
        nft.location = shift_location(ft.location, cut, L)
        new_rec.features.append(nft)

    new_rec.annotations = dict(rec.annotations) if rec.annotations else {}
    if new_rec.annotations.get("topology", "").lower() != "circular":
        new_rec.annotations["topology"] = "circular"

    return new_rec


def rotate_reference_gb_files_to_rrnS(
    ref_dir: Path,
    ref_gb_paths: List[Path],
    out_rot_dir: Path,
    rrns_re: re.Pattern,
) -> Tuple[Path, List[Path]]:
    """
    Write rotated copies of reference GBs into out_rot_dir, preserving folder layout:
      out_rot_dir/<ACC>-<Taxon>/<ACC>.gb

    Returns (rotated_ref_dir, rotated_gb_paths).
    """
    out_rot_dir.mkdir(parents=True, exist_ok=True)

    rotated_paths: List[Path] = []
    for gb in ref_gb_paths:
        try:
            rel = gb.relative_to(ref_dir)
        except Exception:
            rel = Path(gb.name)

        out_gb = out_rot_dir / rel
        out_gb.parent.mkdir(parents=True, exist_ok=True)

        recs = list(SeqIO.parse(str(gb), "genbank"))
        if not recs:
            shutil.copy2(gb, out_gb)
            rotated_paths.append(out_gb)
            continue

        new_recs = []
        any_rot = False
        for rec in recs:
            cut = find_rrnS_cut(rec, rrns_re)
            if cut is None or cut == 0:
                new_recs.append(rec)
            else:
                new_recs.append(rotate_record_to_cut(rec, cut))
                any_rot = True

        SeqIO.write(new_recs, str(out_gb), "genbank")
        rotated_paths.append(out_gb)

        tag = "ROT" if any_rot else "KEEP"
        print(f"[{tag}] REF {gb}  ->  {out_gb}")

    return out_rot_dir, rotated_paths


def _canonicalize_gene_id(s: str) -> str:
    """
    Canonicalize mtDNA gene IDs across different annotation styles so that
    the same gene matches between tracks (e.g., ND1 vs nad1, CYTB vs cob, COI vs cox1).
    This ONLY affects internal link-building IDs, not the plotted feature labels.
    """
    if not s:
        return s

    x = str(s).strip().lower()
    x = re.sub(r"_+", "_", x)

    if x in {"coi", "co1"}:
        x = "cox1"
    if x in {"coii", "co2"}:
        x = "cox2"
    if x in {"coiii", "co3"}:
        x = "cox3"

    m = re.match(r"^(?:nd|nad)([1-6])$", x)
    if m:
        x = f"nad{m.group(1)}"
    if x in {"nd4l", "nad4l"}:
        x = "nad4l"

    if x in {"cytb", "cyt_b", "cytochrome_b"}:
        x = "cob"

    if x in {"atp_6", "atp6"}:
        x = "atp6"
    if x in {"atp_8", "atp8"}:
        x = "atp8"

    if x in {"rrnl", "16s", "16s_rrna", "16s_ribosomal_rna"}:
        x = "rrnl"
    if x in {"rrns", "12s", "12s_rrna", "12s_ribosomal_rna"}:
        x = "rrns"

    if x.startswith("trn"):
        m = re.match(r"^trn([a-z])", x)
        if m:
            x = f"trn{m.group(1)}"

    return x


def _clean_gene_id(s: str) -> str:
    s = s.strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^\w\-\.:()]+", "", s)
    s = re.sub(r"\(.*?\)", "", s)
    s = s.strip("_")
    return _canonicalize_gene_id(s)


def _feature_gene_id(ft) -> Optional[str]:
    q = getattr(ft, "qualifiers", {}) or {}
    for key in ("gene", "locus_tag", "product"):
        if key in q and q[key]:
            val = q[key][0] if isinstance(q[key], list) else str(q[key])
            val = str(val)
            if val.strip():
                return _clean_gene_id(val)
    return None


def extract_gene_hits_from_genbank(gb_path: Path, include_trna_rrna: bool) -> List[GeneHit]:
    hits: List[GeneHit] = []
    for rec in SeqIO.parse(str(gb_path), "genbank"):
        seqid = rec.id
        for ft in rec.features:
            if ft.type not in ("CDS", "tRNA", "rRNA"):
                continue
            if ft.type in ("tRNA", "rRNA") and not include_trna_rrna:
                continue

            gene_id = _feature_gene_id(ft)
            if not gene_id:
                continue

            try:
                start = int(ft.location.start) + 1
                end = int(ft.location.end)
            except Exception:
                continue

            if end <= start:
                continue

            strand_val = ft.location.strand
            strand = 1 if (strand_val is None or strand_val >= 0) else -1

            try:
                seq = str(ft.extract(rec.seq)).upper()
            except Exception:
                continue

            if not seq:
                continue

            hits.append(
                GeneHit(
                    gene_id=gene_id,
                    ftype=ft.type,
                    seqid=seqid,
                    start_1based=start,
                    end_1based=end,
                    strand=strand,
                    seq=seq,
                )
            )
    return hits


def nt_identity_global(seq1: str, seq2: str) -> float:
    if not seq1 or not seq2:
        return 0.0
    alns = pairwise2.align.globalms(seq1, seq2, 1, 0, -1, -0.5, one_alignment_only=True)
    if not alns:
        return 0.0
    a1, a2, _, _, _ = alns[0]
    matches = sum((c1 == c2) and (c1 != "-") for c1, c2 in zip(a1, a2))
    aln_len = max(len(a1), 1)
    return 100.0 * matches / aln_len


def build_best_gene_map(hits: List[GeneHit]) -> Dict[str, GeneHit]:
    m: Dict[str, GeneHit] = {}
    for h in hits:
        L = h.end_1based - h.start_1based + 1
        if h.gene_id not in m:
            m[h.gene_id] = h
        else:
            L0 = m[h.gene_id].end_1based - m[h.gene_id].start_1based + 1
            if L > L0:
                m[h.gene_id] = h
    return m


def add_identity_heatbar(
    fig,
    ax_ref,
    vmin,
    vmax,
    cmap,
    norm,
    title="Identity (%)",
    distance_cm=10.0,
    bar_width_cm=0.55,
    right_margin_cm=0.8,
    title_fontsize=15,
    tick_labelsize=13,
    bar_height_scale=2.5,
):
    """
    Aesthetic vertical heat bar aligned to the main synteny axis and placed
    at a real physical distance from the plot.

    distance_cm controls the gap between the right edge of the synteny plot
    and the left edge of the Identity (%) heat bar.

    bar_height_scale controls heatbar height relative to the main synteny axis.
    Example: 1.0 = same height, 2.5 = 2.5 times taller.
    """
    cm_to_in = 1 / 2.54

    # Store current axes positions in physical inches before resizing the figure.
    old_fig_w, old_fig_h = fig.get_size_inches()
    axes_pos_in = []
    for ax in fig.axes:
        bb = ax.get_position()
        axes_pos_in.append(
            (
                ax,
                bb.x0 * old_fig_w,
                bb.y0 * old_fig_h,
                bb.width * old_fig_w,
                bb.height * old_fig_h,
            )
        )

    ref_bb = ax_ref.get_position()
    ref_x1_in = ref_bb.x1 * old_fig_w
    ref_y0_in = ref_bb.y0 * old_fig_h
    ref_h_in = ref_bb.height * old_fig_h

    pad_in = float(distance_cm) * cm_to_in
    bar_w_in = float(bar_width_cm) * cm_to_in
    right_margin_in = float(right_margin_cm) * cm_to_in

    bar_left_in = ref_x1_in + pad_in
    required_fig_w = bar_left_in + bar_w_in + right_margin_in

    # If needed, enlarge the canvas to keep exactly the requested physical gap.
    new_fig_w = max(old_fig_w, required_fig_w)
    if new_fig_w > old_fig_w:
        fig.set_size_inches(new_fig_w, old_fig_h, forward=True)

        # Keep the synteny plot and existing axes at their original physical size.
        for ax, x0_in, y0_in, w_in, h_in in axes_pos_in:
            ax.set_position([x0_in / new_fig_w, y0_in / old_fig_h, w_in / new_fig_w, h_in / old_fig_h])

    # Make the heatbar taller than the synteny axis, keeping it centered.
    # The coordinates passed to add_axes are still normalized to the figure.
    height_in = ref_h_in * float(bar_height_scale)
    bottom_in = ref_y0_in - (height_in - ref_h_in) / 2

    # Avoid clipping outside the canvas.
    min_margin_in = 0.15
    if bottom_in < min_margin_in:
        bottom_in = min_margin_in
    if bottom_in + height_in > old_fig_h - min_margin_in:
        height_in = max(0.1, old_fig_h - min_margin_in - bottom_in)

    bottom = bottom_in / old_fig_h
    height = height_in / old_fig_h
    left = bar_left_in / new_fig_w
    bar_w = bar_w_in / new_fig_w

    ax_bar = fig.add_axes([left, bottom, bar_w, height])

    grad = np.linspace(vmin, vmax, 512).reshape(-1, 1)
    ax_bar.imshow(
        grad,
        aspect="auto",
        origin="lower",
        cmap=cmap,
        norm=norm,
        extent=(0, 1, vmin, vmax),
    )

    ax_bar.set_xticks([])
    ax_bar.set_xlim(0, 1)

    if (vmax - vmin) >= 15:
        ticks = [vmin, (vmin + vmax) / 2, vmax]
    else:
        ticks = [vmin, vmax]

    ax_bar.set_yticks(ticks)
    ax_bar.set_yticklabels([f"{int(t)}" if float(t).is_integer() else f"{t:.1f}" for t in ticks])

    ax_bar.yaxis.tick_right()
    ax_bar.tick_params(axis="y", length=4.0, width=1.0, labelsize=tick_labelsize, pad=6)

    ax_bar.set_title(title, fontsize=title_fontsize, fontweight="bold", pad=10)

    for spine in ax_bar.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
        spine.set_color("0.25")

    ax_bar.set_facecolor("none")


def run_synteny_gene_links(
    gb_paths: List[Path],
    data_dir: Path,
    out_png: Path,
    out_html: Path,
    identity_thr: float,
    length_thr: int,
    barsize: int,
    dpi: int,
    alpha: float,
    show_trna_rrna: bool,
    ref_dir: Optional[Path] = None,
):
    if len(gb_paths) < 2:
        print(f"[SKIP] Menos de 2 GenBank -> no hago synteny: {out_png.name}", file=sys.stderr)
        return

    gb_paths = sorted(
        gb_paths,
        key=lambda p: sort_key_synteny_custom(p, data_dir, ref_dir=ref_dir),
    )

    gbk_list = []
    for p in gb_paths:
        label = track_label_from_gb(p, data_dir, ref_dir=ref_dir)
        gbk_list.append(Genbank(str(p), name=label))

    gene_hits_per_track: List[Dict[str, GeneHit]] = []
    for p in gb_paths:
        hits = extract_gene_hits_from_genbank(p, include_trna_rrna=show_trna_rrna)
        gene_hits_per_track.append(build_best_gene_map(hits))

    gv = GenomeViz(track_align_type="center", fig_track_height=0.65)

    try:
        try:
            gv.set_scale_bar()
        except TypeError:
            gv.set_scale_bar(barsize=barsize)
    except AttributeError:
        try:
            gv.set_scale_xticks(ymargin=0.6)
        except Exception:
            pass

    for gbk in gbk_list:
        try:
            track = gv.add_feature_track(
                gbk.name,
                gbk.get_seqid2size(),
                labelsize=12,
                align_label=False,
            )
        except TypeError:
            track = gv.add_feature_track(gbk.name, gbk.get_seqid2size())

        for seqid, feats in gbk.get_seqid2features("CDS").items():
            seg = track.get_segment(seqid)
            try:
                seg.add_features(
                    feats,
                    plotstyle="bigarrow",
                    fc=CDS_COLOR,
                    lw=0.8,
                    text_kws=dict(rotation=0, vpos="top", hpos="center", size=9),
                )
            except TypeError:
                for ft in feats:
                    start = int(ft.location.start) + 1
                    end = int(ft.location.end)
                    strand = 1 if (ft.location.strand or 1) >= 0 else -1
                    seg.add_feature(start, end, strand, plotstyle="bigarrow", color=CDS_COLOR, lw=0.8)

        if show_trna_rrna:
            for seqid, feats in gbk.get_seqid2features("tRNA").items():
                seg = track.get_segment(seqid)
                try:
                    seg.add_features(feats, plotstyle="bigarrow", fc=TRNA_COLOR, lw=0.7)
                except TypeError:
                    for ft in feats:
                        start = int(ft.location.start) + 1
                        end = int(ft.location.end)
                        strand = 1 if (ft.location.strand or 1) >= 0 else -1
                        seg.add_feature(start, end, strand, plotstyle="bigarrow", color=TRNA_COLOR, lw=0.7)

            for seqid, feats in gbk.get_seqid2features("rRNA").items():
                seg = track.get_segment(seqid)
                try:
                    seg.add_features(feats, plotstyle="bigarrow", fc=RRNA_COLOR, lw=0.7)
                except TypeError:
                    for ft in feats:
                        start = int(ft.location.start) + 1
                        end = int(ft.location.end)
                        strand = 1 if (ft.location.strand or 1) >= 0 else -1
                        seg.add_feature(start, end, strand, plotstyle="bigarrow", color=RRNA_COLOR, lw=0.7)

    norm = Normalize(vmin=float(identity_thr), vmax=100.0, clip=True)
    cmap = matplotlib.cm.Greys

    any_links = False
    all_identities: List[float] = []

    for i in range(len(gbk_list) - 1):
        left_gbk = gbk_list[i]
        right_gbk = gbk_list[i + 1]

        left_map = gene_hits_per_track[i]
        right_map = gene_hits_per_track[i + 1]

        shared = sorted(set(left_map.keys()) & set(right_map.keys()))
        if not shared:
            continue

        for gene_id in shared:
            Lh = left_map[gene_id]
            Rh = right_map[gene_id]

            min_len = min(len(Lh.seq), len(Rh.seq))
            if min_len < length_thr:
                continue

            ident = nt_identity_global(Lh.seq, Rh.seq)
            if ident < identity_thr:
                continue

            all_identities.append(ident)
            any_links = True

            color = cmap(norm(ident))

            left_link = (left_gbk.name, Lh.seqid, Lh.start_1based, Lh.end_1based)
            right_link = (right_gbk.name, Rh.seqid, Rh.start_1based, Rh.end_1based)

            try:
                gv.add_link(
                    left_link, right_link,
                    color=color,
                    inverted_color="red",
                    curve=True,
                    alpha=alpha,
                )
            except TypeError:
                gv.add_link(left_link, right_link, color=color, alpha=alpha)

    fig = gv.plotfig()

    fig.subplots_adjust(right=0.90)

    if any_links and all_identities:
        try:
            ax0 = fig.axes[0] if fig.axes else None
            if ax0 is not None:
                add_identity_heatbar(
                    fig,
                    ax_ref=ax0,
                    vmin=float(identity_thr),
                    vmax=100.0,
                    cmap=cmap,
                    norm=norm,
                    title="Identity (%)",
                    distance_cm=10.0,
                    bar_width_cm=1.85,
                    title_fontsize=15,
                    tick_labelsize=13,
                    bar_height_scale=5.5,
                )
        except Exception:
            pass

    out_svg = out_png.with_suffix(".svg")

    fig.savefig(out_png, dpi=dpi)
    fig.savefig(out_svg, format="svg")

    try:
        gv.savefig_html(str(out_html))
    except Exception:
        pass

    if any_links:
        print(f"[OK] PNG : {out_png}")
        print(f"[OK] SVG : {out_svg}")
        print(f"[OK] HTML: {out_html}")
    else:
        print(f"[WARN] No links passed the filters. The plot was still saved: {out_png}")
        print(f"[OK] SVG : {out_svg}")


def main():
    a = parse_args()
    root = Path(a.root).expanduser().resolve()
    data_dir = (root / a.data_subdir).resolve()
    outdir = Path(a.outdir).expanduser().resolve() if a.outdir else (root / "results" / "annotation" / "plots" / "synteny").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    by_genus_dir = outdir / "by_genus"
    by_genus_dir.mkdir(parents=True, exist_ok=True)

    genera = list_genera(data_dir)

    ref_dir: Optional[Path] = Path(a.ref_dir).expanduser().resolve()
    ref_gb_paths: List[Path] = []
    if ref_dir.is_dir():
        ref_gb_paths = sorted(ref_dir.rglob(a.ref_pattern))
        if ref_gb_paths:
            print(f"[INFO] References found: N={len(ref_gb_paths)} in {ref_dir}")
        else:
            print(f"[WARN] ref-dir exists but no .gb files matched {a.ref_pattern}: {ref_dir}", file=sys.stderr)
    else:
        print(f"[WARN] ref-dir does not exist: {ref_dir} (continuing without references)", file=sys.stderr)
        ref_dir = None

    if ref_dir and ref_gb_paths and (not a.no_rotate_ref_to_rrnS):
        rrns_re = re.compile(a.rrns_regex, re.IGNORECASE)
        rotated_ref_dir = outdir / "_ref_rotated_rrnS"
        print(f"[INFO] Rotating references to rrnS end + 10 nt -> {rotated_ref_dir}")
        rotated_ref_dir, rotated_ref_gb_paths = rotate_reference_gb_files_to_rrnS(
            ref_dir=ref_dir,
            ref_gb_paths=ref_gb_paths,
            out_rot_dir=rotated_ref_dir,
            rrns_re=rrns_re,
        )
        ref_dir = rotated_ref_dir
        ref_gb_paths = rotated_ref_gb_paths
        print(f"[INFO] Using rotated references: N={len(ref_gb_paths)}")
    elif ref_dir and ref_gb_paths:
        print("[INFO] Reference rotation disabled -> using GenBank files as-is.")

    all_gb: List[Path] = []
    for genus_dir in genera:
        genus = genus_dir.name
        gb_paths = find_gb_in_dir(genus_dir, a.pattern)

        if not gb_paths:
            print(f"[WARN] {genus}: no GenBank files matched {a.pattern}", file=sys.stderr)
            continue

        if ref_dir and ref_gb_paths:
            genus_name = genus_dir.name.lower()
            genus_refs: List[Path] = []
            for rp in ref_gb_paths:
                lbl_plain = track_label_plain_from_gb(rp, data_dir, ref_dir=ref_dir)
                genus_token = lbl_plain.split(" ", 1)[0].strip().lower() if lbl_plain else ""
                if genus_token == genus_name:
                    genus_refs.append(rp)
            if genus_refs:
                gb_paths = gb_paths + genus_refs

        gb_paths = sorted(
            gb_paths,
            key=lambda p: sort_key_synteny_custom(p, data_dir, ref_dir=ref_dir),
        )

        all_gb.extend(gb_paths)

        out_png = by_genus_dir / f"{genus}_{a.prefix_suffix}.png"
        out_html = by_genus_dir / f"{genus}_{a.prefix_suffix}.html"

        print(f"\n[RUN] Genus={genus} | N={len(gb_paths)}")
        run_synteny_gene_links(
            gb_paths=gb_paths,
            data_dir=data_dir,
            out_png=out_png,
            out_html=out_html,
            identity_thr=a.identity_thr,
            length_thr=a.length_thr_nt,
            barsize=a.barsize,
            dpi=a.dpi,
            alpha=a.alpha,
            show_trna_rrna=a.show_trna_rrna,
            ref_dir=ref_dir,
        )

    if ref_dir and ref_gb_paths:
        all_gb.extend(ref_gb_paths)

    all_gb = sorted(
        set(all_gb),
        key=lambda p: sort_key_synteny_custom(p, data_dir, ref_dir=ref_dir),
    )
    if len(all_gb) >= 2:
        out_png = outdir / f"{a.prefix_all}.png"
        out_html = outdir / f"{a.prefix_all}.html"
        print(f"\n[RUN] GLOBAL | N={len(all_gb)}")
        run_synteny_gene_links(
            gb_paths=all_gb,
            data_dir=data_dir,
            out_png=out_png,
            out_html=out_html,
            identity_thr=a.identity_thr,
            length_thr=a.length_thr_nt,
            barsize=a.barsize,
            dpi=a.dpi,
            alpha=a.alpha,
            show_trna_rrna=a.show_trna_rrna,
            ref_dir=ref_dir,
        )
    else:
        print("[WARN] Global: fewer than 2 GenBank files, skipping global synteny.", file=sys.stderr)


if __name__ == "__main__":
    main()
