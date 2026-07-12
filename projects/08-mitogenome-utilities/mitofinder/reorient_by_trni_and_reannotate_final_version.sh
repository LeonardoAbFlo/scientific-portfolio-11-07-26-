#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ============================================================
# Medaka consensus -> orientar por trnI -> correr MitoFinder
# Maneja barcodes invertidos usando la hebra de trnI.
# ============================================================

RESULTS_P1="$HOME/mosquitos/results_p1"

DATA_ROOT="$RESULTS_P1/data"

MEDAKA_ORG="$RESULTS_P1/mitogenomes_medaka"
PREV_MITOFINDER="$RESULTS_P1/medaka_mitofinder"

OUT_TRNI="$RESULTS_P1/mitogenomas_medaka_trnI"
OUT_MITOFINDER_TRNI="$RESULTS_P1/medaka_mitofinder_trnI"

REF_DIR="$RESULTS_P1/ref_gb"
REF_GB="$REF_DIR/culicidae_mt_refseq.gb"
REF_GB_ORIGINAL="/path/to/Mosquito_Native_F1/scripts/mosquitos_genbank/culicidae_mt_refseq.gb"

GENCODE=5
THREADS=32
MAXMEM=16

RUN_MITOFINDER=1

REPORT="$RESULTS_P1/medaka_trnI_rotation_report.tsv"

mkdir -p "$REF_DIR" "$OUT_TRNI" "$OUT_MITOFINDER_TRNI"

if [[ ! -f "$REF_GB" ]]; then
  if [[ -f "$REF_GB_ORIGINAL" ]]; then
    cp "$REF_GB_ORIGINAL" "$REF_GB"
  else
    echo "ERROR: No encontré la referencia GenBank:"
    echo "  $REF_GB"
    echo "  $REF_GB_ORIGINAL"
    exit 1
  fi
fi

if [[ ! -d "$DATA_ROOT" ]]; then
  echo "ERROR: No existe DATA_ROOT: $DATA_ROOT"
  exit 1
fi

if [[ ! -d "$MEDAKA_ORG" ]]; then
  echo "ERROR: No existe MEDAKA_ORG: $MEDAKA_ORG"
  exit 1
fi

if [[ ! -d "$PREV_MITOFINDER" ]]; then
  echo "ERROR: No existe PREV_MITOFINDER: $PREV_MITOFINDER"
  exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate sif_env

python3 - <<'PY'
try:
    import Bio
except Exception:
    raise SystemExit("ERROR: Biopython no está instalado. Usa: conda install -c conda-forge biopython")
PY

cat > "$REPORT" <<'EOF'
sample	genus	input_fasta	genbank_used	trnI_start_1based	trnI_end_1based	trnI_strand	len_input_fasta	len_genbank	sequence_used	operation	output_fasta	status
EOF

echo "============================================================"
echo "Orientando mitogenomas Medaka por trnI"
echo "Maneja barcodes invertidos con reverse-complement"
echo "OUT_TRNI:            $OUT_TRNI"
echo "OUT_MITOFINDER_TRNI: $OUT_MITOFINDER_TRNI"
echo "REPORT:              $REPORT"
echo "============================================================"
echo

processed=0
rotated_ok=0
reverse_ok=0
trni_not_found=0
missing_consensus=0
missing_genbank=0
failed_mitofinder=0

for GENUS_DIR in "$DATA_ROOT"/*/; do
  GENUS="$(basename "$GENUS_DIR")"

  if [[ "$GENUS" == "ref_gb" ]]; then
    continue
  fi

  for SAMPLE_DIR in "$GENUS_DIR"*/; do
    SAMPLE_NAME="$(basename "$SAMPLE_DIR")"

    if [[ "$SAMPLE_NAME" =~ ^(F[0-9]+)_([0-9]+)-(.+)$ ]]; then
      FLOWCELL="${BASH_REMATCH[1]}"
      BARCODE_NUM="${BASH_REMATCH[2]}"
      BARCODE_DIR="barcode${BARCODE_NUM}"
    else
      echo "ADVERTENCIA: No pude interpretar: $SAMPLE_NAME"
      continue
    fi

    SAMPLE_SAFE="$(echo "$SAMPLE_NAME" | sed -E 's/[^A-Za-z0-9_]+/_/g; s/_+$//')"

    CONSENSUS="$MEDAKA_ORG/$GENUS/$SAMPLE_NAME/consensus.fasta"
    PREV_MF_DIR="$PREV_MITOFINDER/$GENUS/$SAMPLE_NAME"

    TRNI_SAMPLE_DIR="$OUT_TRNI/$GENUS/$SAMPLE_NAME"
    TRNI_FASTA="$TRNI_SAMPLE_DIR/${SAMPLE_SAFE}_trnI.fasta"

    mkdir -p "$TRNI_SAMPLE_DIR"

    echo "------------------------------------------------------------"
    echo "Muestra:    $SAMPLE_NAME"
    echo "Género:     $GENUS"
    echo "Flowcell:   $FLOWCELL"
    echo "Barcode:    $BARCODE_DIR"
    echo "Consensus:  $CONSENSUS"
    echo "Salida:     $TRNI_FASTA"

    if [[ ! -s "$CONSENSUS" ]]; then
      echo "  ! No existe consensus.fasta"
      echo -e "${SAMPLE_NAME}\t${GENUS}\t${CONSENSUS}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t${TRNI_FASTA}\tMISSING_CONSENSUS" >> "$REPORT"
      ((missing_consensus+=1))
      echo
      continue
    fi

    if [[ ! -d "$PREV_MF_DIR" ]]; then
      echo "  ! No existe MitoFinder previo: $PREV_MF_DIR"
      cp "$CONSENSUS" "$TRNI_FASTA"
      echo -e "${SAMPLE_NAME}\t${GENUS}\t${CONSENSUS}\tNA\tNA\tNA\tNA\tNA\tNA\tinput_fasta_unrotated\tcopy_unrotated\t${TRNI_FASTA}\tMISSING_PREVIOUS_MITOFINDER" >> "$REPORT"
      ((missing_genbank+=1))
    else
      mapfile -t GB_FILES < <(
        find "$PREV_MF_DIR" -type f \
          \( -name "*.gb" -o -name "*.gbk" -o -name "*.gbff" \) \
          2>/dev/null | sort
      )

      if ((${#GB_FILES[@]} == 0)); then
        echo "  ! No encontré GenBank previo"
        cp "$CONSENSUS" "$TRNI_FASTA"
        echo -e "${SAMPLE_NAME}\t${GENUS}\t${CONSENSUS}\tNA\tNA\tNA\tNA\tNA\tNA\tinput_fasta_unrotated\tcopy_unrotated\t${TRNI_FASTA}\tMISSING_GENBANK" >> "$REPORT"
        ((missing_genbank+=1))
      else
        echo "  → Buscando trnI y orientación..."

        if python3 - "$CONSENSUS" "$TRNI_FASTA" "$SAMPLE_NAME" "$GENUS" "${GB_FILES[@]}" >> "$REPORT" <<'PY'
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import re

input_fasta = Path(sys.argv[1])
output_fasta = Path(sys.argv[2])
sample = sys.argv[3]
genus = sys.argv[4]
gb_files = [Path(x) for x in sys.argv[5:]]

def feature_text(feat):
    parts = []
    for key, values in feat.qualifiers.items():
        if isinstance(values, list):
            parts.extend(str(v) for v in values)
        else:
            parts.append(str(values))
    return " ".join(parts)

def is_trnI(feat):
    ftype = feat.type.lower()
    if ftype not in {"trna", "gene"}:
        return False

    txt = feature_text(feat)

    patterns = [
        r"\btrnI\b",
        r"\btrn-I\b",
        r"\btrn_I\b",
        r"trna[-_ ]?ile",
        r"tRNA[-_ ]?Ile",
        r"transfer RNA[-_ ]?Ile",
        r"isoleucine",
    ]

    return any(re.search(p, txt, flags=re.IGNORECASE) for p in patterns)

def location_bounds_and_strand(feat):
    loc = feat.location

    if hasattr(loc, "parts"):
        starts = [int(p.start) for p in loc.parts]
        ends = [int(p.end) for p in loc.parts]
        strands = [p.strand for p in loc.parts if p.strand is not None]
    else:
        starts = [int(loc.start)]
        ends = [int(loc.end)]
        strands = [loc.strand] if loc.strand is not None else []

    start = min(starts)
    end = max(ends)

    if -1 in strands and 1 not in strands:
        strand = -1
    elif 1 in strands and -1 not in strands:
        strand = 1
    else:
        strand = getattr(loc, "strand", None)
        if strand is None:
            strand = 1

    return start, end, strand

fasta_records = list(SeqIO.parse(str(input_fasta), "fasta"))
if not fasta_records:
    raise SystemExit("NO_FASTA_RECORDS")

input_rec = fasta_records[0]
input_seq = input_rec.seq

selected = None

# Preferir features tRNA sobre gene
candidate_features = []

for gb in gb_files:
    try:
        for rec in SeqIO.parse(str(gb), "genbank"):
            for feat in rec.features:
                if is_trnI(feat):
                    candidate_features.append((gb, rec, feat))
    except Exception:
        continue

if not candidate_features:
    raise SystemExit("TRNI_NOT_FOUND")

candidate_features.sort(key=lambda x: 0 if x[2].type.lower() == "trna" else 1)
gb_used, gb_rec, feat = candidate_features[0]

start, end, strand = location_bounds_and_strand(feat)

gb_seq = gb_rec.seq

if len(input_seq) == len(gb_seq):
    seq_to_use = input_seq
    sequence_used = "input_fasta"
else:
    seq_to_use = gb_seq
    sequence_used = "genbank_sequence"

L = len(seq_to_use)

if start < 0 or end > L or start >= end:
    raise SystemExit(f"BAD_TRNI_COORDINATES:{start}:{end}:len={L}")

if strand == -1:
    # Si trnI está en hebra negativa, primero hacemos reverse-complement.
    # Coordenadas BioPython son 0-based semiabiertas [start, end).
    # En la secuencia reverse-complement, el inicio equivalente es L - end.
    rc_seq = seq_to_use.reverse_complement()
    rc_start = L - end
    rotated_seq = rc_seq[rc_start:] + rc_seq[:rc_start]
    operation = "reverse_complement_then_rotate"
else:
    rotated_seq = seq_to_use[start:] + seq_to_use[:start]
    operation = "rotate_only"

out_rec = SeqRecord(
    rotated_seq,
    id=f"{sample}_trnI",
    name=f"{sample}_trnI",
    description=(
        f"oriented_to_start_at_trnI "
        f"trnI_start_1based={start+1} "
        f"trnI_end_1based={end} "
        f"strand={strand} "
        f"operation={operation} "
        f"source={sequence_used}"
    )
)

output_fasta.parent.mkdir(parents=True, exist_ok=True)
SeqIO.write(out_rec, str(output_fasta), "fasta")

print(
    f"{sample}\t{genus}\t{input_fasta}\t{gb_used}\t"
    f"{start+1}\t{end}\t{strand}\t"
    f"{len(input_seq)}\t{len(gb_seq)}\t"
    f"{sequence_used}\t{operation}\t{output_fasta}\tROTATED_OK"
)
PY
        then
          echo "  ✓ FASTA orientado por trnI"

          if tail -n 1 "$REPORT" | grep -q "reverse_complement_then_rotate"; then
            echo "  ✓ Barcode invertido corregido con reverse-complement"
            ((reverse_ok+=1))
          fi

          ((rotated_ok+=1))
        else
          echo "  ! No se pudo ubicar trnI. Se copiará sin rotar."
          cp "$CONSENSUS" "$TRNI_FASTA"
          echo -e "${SAMPLE_NAME}\t${GENUS}\t${CONSENSUS}\tNA\tNA\tNA\tNA\tNA\tNA\tinput_fasta_unrotated\tcopy_unrotated\t${TRNI_FASTA}\tTRNI_NOT_FOUND_OR_FAILED" >> "$REPORT"
          ((trni_not_found+=1))
        fi
      fi
    fi

    if [[ ! -s "$TRNI_FASTA" ]]; then
      echo "  ! No se generó FASTA válido"
      echo
      continue
    fi

    if [[ "$RUN_MITOFINDER" == "1" ]]; then
      MF_TRNI_OUTDIR="$OUT_MITOFINDER_TRNI/$GENUS/$SAMPLE_NAME"
      mkdir -p "$MF_TRNI_OUTDIR"

      SAMPLE_TRNI_ID="${SAMPLE_SAFE}_trnI"

      echo "  → Ejecutando MitoFinder sobre FASTA orientado"
      echo "    Input:  $TRNI_FASTA"
      echo "    Output: $MF_TRNI_OUTDIR"

      if ! (
        cd "$MF_TRNI_OUTDIR"

        mitofinder \
          --seqid "$SAMPLE_TRNI_ID" \
          --assembly "$TRNI_FASTA" \
          --refseq "$REF_GB" \
          --organism "$GENCODE" \
          --processors "$THREADS" \
          --max-memory "$MAXMEM" \
          --override \
          --blast-size 30 \
          --min-contig-size 1000
      ); then
        echo "  ! ERROR: MitoFinder falló para $SAMPLE_NAME"
        ((failed_mitofinder+=1))
        echo
        continue
      fi

      COI_FASTA="$MF_TRNI_OUTDIR/COI_${SAMPLE_TRNI_ID}.fna"

      mapfile -t CDS_FILES < <(
        find "$MF_TRNI_OUTDIR" -type f -path "*/CDS/*.fna" 2>/dev/null | sort
      )

      if ((${#CDS_FILES[@]} > 0)); then
        awk '
          BEGIN { IGNORECASE=1 }
          /^>/ {
            keep = ($0 ~ /(COI|COX1|CO1|cytochrome c oxidase subunit I|cytochrome oxidase subunit I)/)
          }
          keep { print }
        ' "${CDS_FILES[@]}" > "${COI_FASTA}.tmp"

        if [[ -s "${COI_FASTA}.tmp" ]]; then
          mv "${COI_FASTA}.tmp" "$COI_FASTA"
          echo "  ✓ COI/COX1 guardado en: $COI_FASTA"
        else
          rm -f "${COI_FASTA}.tmp"
          echo "  ! No se encontró COI/COX1 en la corrida trnI"
        fi
      else
        echo "  ! No se encontraron CDS/*.fna en la corrida trnI"
      fi
    fi

    ((processed+=1))
    echo

  done
done

echo "============================================================"
echo "PROCESO COMPLETADO"
echo "Muestras procesadas:         $processed"
echo "Reordenadas por trnI:        $rotated_ok"
echo "Invertidas corregidas:       $reverse_ok"
echo "trnI no encontrado/fallido:  $trni_not_found"
echo "Consensus faltantes:         $missing_consensus"
echo "GenBank previos faltantes:   $missing_genbank"
echo "MitoFinder fallidos:         $failed_mitofinder"
echo
echo "FASTA orientados en:"
echo "  $OUT_TRNI"
echo
echo "MitoFinder trnI en:"
echo "  $OUT_MITOFINDER_TRNI"
echo
echo "Reporte:"
echo "  $REPORT"
echo "============================================================"