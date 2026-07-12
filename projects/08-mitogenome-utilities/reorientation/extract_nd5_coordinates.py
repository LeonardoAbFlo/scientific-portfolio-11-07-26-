#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def is_nd5(gene: str) -> bool:
    g = gene.strip().lower()
    return g.startswith("nad5") or g.startswith("nd5")


def parse_mitos2_faa_header(line: str) -> Optional[Tuple[str, int, int, str, str]]:
    """
    Parsea headers tipo:
      >barcode95; 6362-6922; -; nad5_1
    Retorna (barcode, start, end, strand, gene)
    """
    if not line.startswith(">"):
        return None
    s = line[1:].strip()

    # Split por ';' (MITOS2 usa ese formato)
    parts = [p.strip() for p in s.split(";")]
    if len(parts) < 4:
        return None

    barcode = parts[0]
    coord = parts[1].replace(" ", "")
    strand = parts[2]
    gene = parts[3].split()[0]  # por si hay texto extra

    m = re.match(r"^(\d+)-(\d+)$", coord)
    if not m:
        return None

    start = int(m.group(1))
    end = int(m.group(2))

    if strand not in {"+", "-"}:
        return None

    return barcode, start, end, strand, gene


def nd5_sort_key(start: int, end: int, strand: str) -> Tuple[int, int]:
    """
    Orden lógico de fragmentos ND5:
      - strand '+' : por start ascendente
      - strand '-' : por end descendente (5' suele estar en coordenada alta)
    """
    if strand == "-":
        return (-end, -start)
    return (start, end)


def iter_result_faa(root: Path):
    """
    Busca result.faa en patrón root/<genus>/<sample>/result.faa
    Ignora _logs.
    """
    for genus_dir in sorted([p for p in root.iterdir() if p.is_dir()]):
        if genus_dir.name == "_logs":
            continue
        for sample_dir in sorted([p for p in genus_dir.iterdir() if p.is_dir()]):
            faa = sample_dir / "result.faa"
            if faa.is_file():
                yield genus_dir.name, sample_dir.name, faa


def main():
    ap = argparse.ArgumentParser(
        description="Extrae coordenadas ND5 (MITOS2 result.faa) por muestra y genera TSV/JSON."
    )
    ap.add_argument(
        "--root",
        required=True,
        nargs="+",
        help="Ruta(s) a mitos2_out (ej: ~/mosquitos/results_p1/mitos2_out).",
    )
    ap.add_argument(
        "--out-tsv",
        required=True,
        help="Salida TSV con ND5 por muestra.",
    )
    ap.add_argument(
        "--out-json",
        default=None,
        help="(Opcional) Salida JSON con ND5 por muestra.",
    )
    args = ap.parse_args()

    roots = [Path(r).expanduser().resolve() for r in args.root]
    out_tsv = Path(args.out_tsv).expanduser().resolve()
    out_json = Path(args.out_json).expanduser().resolve() if args.out_json else None

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    if out_json:
        out_json.parent.mkdir(parents=True, exist_ok=True)

    rows: List[str] = []
    rows.append("\t".join([
        "sample_id",
        "genus_dir",
        "barcode",
        "strand",
        "parts",
        "coords_parts_ordered",
        "win_start",
        "win_end",
        "source_result_faa",
    ]))

    data: Dict[str, dict] = {}
    n = 0

    for root in roots:
        if not root.is_dir():
            print(f"[WARN] root inválido: {root}", file=sys.stderr)
            continue

        for genus_dir, sample_id, faa_path in iter_result_faa(root):
            nd5_parts = []
            barcodes = set()
            strands = set()

            try:
                with faa_path.open("r", encoding="utf-8", errors="replace") as f:
                    for line in f:
                        if not line.startswith(">"):
                            continue
                        parsed = parse_mitos2_faa_header(line)
                        if not parsed:
                            continue
                        barcode, start, end, strand, gene = parsed
                        if not is_nd5(gene):
                            continue

                        nd5_parts.append({
                            "barcode": barcode,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "gene": gene,
                        })
                        barcodes.add(barcode)
                        strands.add(strand)
            except Exception as e:
                print(f"[WARN] no pude leer {faa_path}: {e}", file=sys.stderr)
                continue

            if not nd5_parts:
                # no ND5 para esa muestra
                continue

            # Si hay inconsistencias, lo dejamos señalado pero igual exportamos
            barcode_show = sorted(barcodes)[0] if barcodes else "NA"
            strand_show = sorted(strands)[0] if strands else "NA"
            if len(barcodes) > 1:
                print(f"[WARN] múltiples barcodes en {sample_id}: {sorted(barcodes)}", file=sys.stderr)
            if len(strands) > 1:
                print(f"[WARN] múltiples strands en {sample_id}: {sorted(strands)}", file=sys.stderr)

            # Ordenar fragmentos en orden "biológico" según hebra
            nd5_parts.sort(key=lambda p: nd5_sort_key(p["start"], p["end"], p["strand"]))

            # Ventana ND5: mínimo start y máximo end (coordenadas NT en el contig)
            win_start = min(p["start"] for p in nd5_parts)
            win_end = max(p["end"] for p in nd5_parts)

            coords_parts = ",".join([f'{p["start"]}-{p["end"]}:{p["gene"]}' for p in nd5_parts])

            rows.append("\t".join([
                sample_id,
                genus_dir,
                barcode_show,
                strand_show,
                str(len(nd5_parts)),
                coords_parts,
                str(win_start),
                str(win_end),
                str(faa_path),
            ]))

            data[sample_id] = {
                "genus_dir": genus_dir,
                "barcode": barcode_show,
                "strand": strand_show,
                "parts": nd5_parts,
                "win_start": win_start,
                "win_end": win_end,
                "source_result_faa": str(faa_path),
            }

            n += 1

    out_tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"[OK] TSV escrito: {out_tsv} | muestras con ND5: {n}")

    if out_json:
        out_json.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[OK] JSON escrito: {out_json}")


if __name__ == "__main__":
    main()


"""
Usage:
./extract_nd5_coordinates.py \
  --root ~/mosquitos/results_p1/mitos2_out \
  --out-tsv ~/mosquitos/results_p1/ND5_coords.tsv \
  --out-json ~/mosquitos/results_p1/ND5_coords.json

"""
