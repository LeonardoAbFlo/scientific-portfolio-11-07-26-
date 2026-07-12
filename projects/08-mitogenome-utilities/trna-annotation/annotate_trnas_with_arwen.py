#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import subprocess
from Bio import SeqIO
import pandas as pd

# === CONFIGURACIÓN ===
PROJECT_ROOT = "/path/to/mosquitos/F1"
MITOF_DIR    = os.path.join(PROJECT_ROOT, "mitofinder_trnI")

# Carpeta donde quedará todo:
OUTDIR = os.path.join(PROJECT_ROOT, "output", "tRNA_sec_struc")
os.makedirs(OUTDIR, exist_ok=True)

ARWEN_BIN = "arwen"   # asumiendo que está en el PATH / conda env activo
# =====================

def extract_trnas_from_gb(gb_path, fasta_out, tsv_out):
    """Extrae tRNAs de un GenBank y escribe FASTA + TSV."""
    rec = SeqIO.read(gb_path, "genbank")

    rows = []
    with open(fasta_out, "w") as outfa:
        for f in rec.features:
            if f.type == "tRNA":
                seq = f.location.extract(rec.seq)
                start_1based = int(f.location.start) + 1
                end_1based   = int(f.location.end)   # end exclusivo en Biopython
                strand = "+" if f.location.strand == 1 else "-"
                name = (f.qualifiers.get("gene") or
                        f.qualifiers.get("product") or
                        ["tRNA"])[0]
                anticodon = (f.qualifiers.get("anticodon") or [""])[0]

                ident = f"{name}|{start_1based}-{end_1based}({strand})"
                outfa.write(f">{ident}\n{str(seq)}\n")

                rows.append([
                    name,
                    start_1based,
                    end_1based,
                    strand,
                    anticodon,
                    len(seq)
                ])

    df = pd.DataFrame(
        rows,
        columns=["tRNA", "start", "end", "strand", "anticodon", "length"]
    )
    df.to_csv(tsv_out, sep="\t", index=False)


def main():
    # Busca todos los GenBank tipo:
    # mitofinder_trnI/barcode81_trnI/barcode81_trnI_MitoFinder_mitfi_Final_Results/barcode81_trnI_mtDNA_contig.gb
    pattern = os.path.join(
        MITOF_DIR,
        "barcode*_trnI",
        "*_MitoFinder_mitfi_Final_Results",
        "*_mtDNA_contig.gb"
    )
    gb_files = sorted(glob.glob(pattern))

    if not gb_files:
        print(f"[ERROR] No encontré GenBank con patrón: {pattern}")
        return

    print(f"[INFO] Encontrados {len(gb_files)} GenBank.")

    for gb_path in gb_files:
        gb_base = os.path.basename(gb_path)  # ej: barcode81_trnI_mtDNA_contig.gb
        # Nos quedamos con "barcode81" como nombre de barcode
        barcode = gb_base.split("_")[0]

        print(f"\n[INFO] Procesando {barcode}")
        print(f"       GenBank: {gb_path}")

        # Prefijo común de salida para este barcode
        prefix = os.path.join(OUTDIR, barcode)

        fasta_out = prefix + "_tRNAs.fasta"
        tsv_out   = prefix + "_tRNAs_coordinates.tsv"
        arwen_out = prefix + ".out"   # para que el nombre sea, p.ej., barcode81.out

        # 1) Extraer tRNAs
        extract_trnas_from_gb(gb_path, fasta_out, tsv_out)
        print(f"[OK] tRNAs FASTA: {fasta_out}")
        print(f"[OK] Coordenadas: {tsv_out}")

        # 2) Ejecutar ARWEN
        #    Equivalente a: arwen -sequence tRNAs.fasta -mt -o barcode81.out
        cmd = [
            ARWEN_BIN,
            "-sequence", fasta_out,
            "-mt",
            "-o", arwen_out
        ]
        print(f"[CMD] {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True)
            print(f"[OK] ARWEN output (prefijo): {arwen_out}")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Falló ARWEN para {barcode}: {e}")

    print(f"\n[FIN] Todos los resultados en: {OUTDIR}")


if __name__ == "__main__":
    main()
