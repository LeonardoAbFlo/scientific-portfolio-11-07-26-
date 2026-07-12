#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ============================================================
# Medaka consensus -> orientar 10 nt después del final de rrnS -> correr MitoFinder
# Fuerza el orden:
# región control -> tRNA-Ile -> tRNA-Gln -> tRNA-Met -> demás genes -> rrnS
# Maneja barcodes invertidos usando reverse-complement cuando sea necesario.
# ============================================================

RESULTS_P1="$HOME/mosquitos/results_p1"

DATA_ROOT="$RESULTS_P1/data"

MEDAKA_ORG="$RESULTS_P1/mitogenomes_medaka"
PREV_MITOFINDER="$RESULTS_P1/medaka_mitofinder"

OUT_RRNS="$RESULTS_P1/mitogenomas_medaka_rrnS"
OUT_MITOFINDER_RRNS="$RESULTS_P1/medaka_mitofinder_rrnS"

REF_DIR="$RESULTS_P1/ref_gb"
REF_GB="$REF_DIR/culicidae_mt_refseq.gb"
REF_GB_ORIGINAL="/path/to/Mosquito_Native_F1/scripts/mosquitos_genbank/culicidae_mt_refseq.gb"

GENCODE=5
THREADS=32
MAXMEM=16

RUN_MITOFINDER=1

REPORT="$RESULTS_P1/medaka_rrnS_rotation_report.tsv"

mkdir -p "$REF_DIR" "$OUT_RRNS" "$OUT_MITOFINDER_RRNS"

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
sample	genus	input_fasta	genbank_used	rrnS_start_1based	rrnS_end_1based	rrnS_strand	cut_position_1based	len_input_fasta	len_genbank	sequence_used	operation	output_fasta	status
EOF

echo "============================================================"
echo "Orientando mitogenomas Medaka 10 nt después del final de rrnS"
echo "Forzando orden: región control -> trnI -> trnQ -> trnM -> ... -> rrnS"
echo "Maneja barcodes invertidos con reverse-complement"
echo "OUT_RRNS:             $OUT_RRNS"
echo "OUT_MITOFINDER_RRNS:  $OUT_MITOFINDER_RRNS"
echo "REPORT:               $REPORT"
echo "============================================================"
echo

processed=0
rotated_ok=0
reverse_ok=0
rrns_not_found=0
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

    RRNS_SAMPLE_DIR="$OUT_RRNS/$GENUS/$SAMPLE_NAME"
    RRNS_FASTA="$RRNS_SAMPLE_DIR/${SAMPLE_SAFE}_rrnS.fasta"

    mkdir -p "$RRNS_SAMPLE_DIR"

    echo "------------------------------------------------------------"
    echo "Muestra:    $SAMPLE_NAME"
    echo "Género:     $GENUS"
    echo "Flowcell:   $FLOWCELL"
    echo "Barcode:    $BARCODE_DIR"
    echo "Consensus:  $CONSENSUS"
    echo "Salida:     $RRNS_FASTA"

    if [[ ! -s "$CONSENSUS" ]]; then
      echo "  ! No existe consensus.fasta"
      echo -e "${SAMPLE_NAME}\t${GENUS}\t${CONSENSUS}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t${RRNS_FASTA}\tMISSING_CONSENSUS" >> "$REPORT"
      ((missing_consensus+=1))
      echo
      continue
    fi

    if [[ ! -d "$PREV_MF_DIR" ]]; then
      echo "  ! No existe MitoFinder previo: $PREV_MF_DIR"
      cp "$CONSENSUS" "$RRNS_FASTA"
      echo -e "${SAMPLE_NAME}\t${GENUS}\t${CONSENSUS}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tinput_fasta_unrotated\tcopy_unrotated\t${RRNS_FASTA}\tMISSING_PREVIOUS_MITOFINDER" >> "$REPORT"
      ((missing_genbank+=1))
    else
      mapfile -t GB_FILES < <(
        find "$PREV_MF_DIR" -type f \
          \( -name "*.gb" -o -name "*.gbk" -o -name "*.gbff" \) \
          2>/dev/null | sort
      )

      if ((${#GB_FILES[@]} == 0)); then
        echo "  ! No encontré GenBank previo"
        cp "$CONSENSUS" "$RRNS_FASTA"
        echo -e "${SAMPLE_NAME}\t${GENUS}\t${CONSENSUS}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tinput_fasta_unrotated\tcopy_unrotated\t${RRNS_FASTA}\tMISSING_GENBANK" >> "$REPORT"
        ((missing_genbank+=1))
      else
        echo "  → Buscando rrnS y corrigiendo orientación canónica..."

        if python3 - "$CONSENSUS" "$RRNS_FASTA" "$SAMPLE_NAME" "$GENUS" "${GB_FILES[@]}" >> "$REPORT" <<'PY'
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

OFFSET_AFTER_RRNS = 10

def feature_text(feat):
    parts = []
    for key, values in feat.qualifiers.items():
        if isinstance(values, list):
            parts.extend(str(v) for v in values)
        else:
            parts.append(str(values))
    return " ".join(parts)

def match_any(txt, patterns):
    return any(re.search(p, txt, flags=re.IGNORECASE) for p in patterns)

def is_rrnS(feat):
    ftype = feat.type.lower()
    if ftype not in {"rrna", "gene"}:
        return False

    txt = feature_text(feat)

    patterns = [
        r"\brrnS\b",
        r"\brrn-S\b",
        r"\brrn_S\b",
        r"\b12S\b",
        r"12S[-_ ]?rRNA",
        r"small[-_ ]?subunit[-_ ]?ribosomal[-_ ]?RNA",
        r"small[-_ ]?subunit[-_ ]?rRNA",
        r"s[-_ ]?rRNA",
        r"srRNA",
    ]

    return match_any(txt, patterns)

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

    return match_any(txt, patterns)

def is_trnQ(feat):
    ftype = feat.type.lower()
    if ftype not in {"trna", "gene"}:
        return False

    txt = feature_text(feat)

    patterns = [
        r"\btrnQ\b",
        r"\btrn-Q\b",
        r"\btrn_Q\b",
        r"trna[-_ ]?gln",
        r"tRNA[-_ ]?Gln",
        r"transfer RNA[-_ ]?Gln",
        r"glutamine",
    ]

    return match_any(txt, patterns)

def is_trnM(feat):
    ftype = feat.type.lower()
    if ftype not in {"trna", "gene"}:
        return False

    txt = feature_text(feat)

    patterns = [
        r"\btrnM\b",
        r"\btrn-M\b",
        r"\btrn_M\b",
        r"trna[-_ ]?met",
        r"tRNA[-_ ]?Met",
        r"transfer RNA[-_ ]?Met",
        r"methionine",
    ]

    return match_any(txt, patterns)

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

def pick_best_feature(features, matcher, preferred_type=None):
    candidates = [f for f in features if matcher(f)]
    if not candidates:
        return None

    if preferred_type is not None:
        candidates.sort(key=lambda f: 0 if f.type.lower() == preferred_type else 1)

    return candidates[0]

def rc_start_of_feature(feat, L):
    start, end, strand = location_bounds_and_strand(feat)
    return L - end

def forward_start_of_feature(feat):
    start, end, strand = location_bounds_and_strand(feat)
    return start

def relative_position(pos, cut_pos, L):
    return (pos - cut_pos) % L

def canonical_order_after_cut(cut_pos, marker_features, L, orientation):
    """
    Retorna True si después del corte aparecen:
    trnI -> trnQ -> trnM
    en esa orientación.
    """

    trnI = marker_features.get("trnI")
    trnQ = marker_features.get("trnQ")
    trnM = marker_features.get("trnM")

    if trnI is None or trnQ is None or trnM is None:
        return False, None

    if orientation == "forward":
        pos_I = forward_start_of_feature(trnI)
        pos_Q = forward_start_of_feature(trnQ)
        pos_M = forward_start_of_feature(trnM)
    elif orientation == "reverse":
        pos_I = rc_start_of_feature(trnI, L)
        pos_Q = rc_start_of_feature(trnQ, L)
        pos_M = rc_start_of_feature(trnM, L)
    else:
        raise ValueError("orientation must be forward or reverse")

    rel_I = relative_position(pos_I, cut_pos, L)
    rel_Q = relative_position(pos_Q, cut_pos, L)
    rel_M = relative_position(pos_M, cut_pos, L)

    order = {
        "trnI": rel_I,
        "trnQ": rel_Q,
        "trnM": rel_M,
    }

    return rel_I < rel_Q < rel_M, order

fasta_records = list(SeqIO.parse(str(input_fasta), "fasta"))
if not fasta_records:
    raise SystemExit("NO_FASTA_RECORDS")

input_rec = fasta_records[0]
input_seq = input_rec.seq

rrns_candidates = []

for gb in gb_files:
    try:
        for rec in SeqIO.parse(str(gb), "genbank"):
            for feat in rec.features:
                if is_rrnS(feat):
                    rrns_candidates.append((gb, rec, feat))
    except Exception:
        continue

if not rrns_candidates:
    raise SystemExit("RRNS_NOT_FOUND")

rrns_candidates.sort(key=lambda x: 0 if x[2].type.lower() == "rrna" else 1)
gb_used, gb_rec, rrns_feat = rrns_candidates[0]

gb_seq = gb_rec.seq

if len(input_seq) == len(gb_seq):
    seq_to_use = input_seq
    sequence_used = "input_fasta"
else:
    seq_to_use = gb_seq
    sequence_used = "genbank_sequence"

L = len(seq_to_use)

start, end, strand = location_bounds_and_strand(rrns_feat)

if start < 0 or end > L or start >= end:
    raise SystemExit(f"BAD_RRNS_COORDINATES:{start}:{end}:len={L}")

all_features = list(gb_rec.features)

marker_features = {
    "trnI": pick_best_feature(all_features, is_trnI, preferred_type="trna"),
    "trnQ": pick_best_feature(all_features, is_trnQ, preferred_type="trna"),
    "trnM": pick_best_feature(all_features, is_trnM, preferred_type="trna"),
}

# ------------------------------------------------------------
# Orientación forward:
# cortar 10 nt después del final de rrnS.
# ------------------------------------------------------------

cut_forward = (end + OFFSET_AFTER_RRNS) % L

forward_ok, forward_order = canonical_order_after_cut(
    cut_forward,
    marker_features,
    L,
    orientation="forward"
)

# ------------------------------------------------------------
# Orientación reverse-complement:
# rrnS en RC tendrá coordenadas:
# start_rc = L - end
# end_rc   = L - start
# cortar 10 nt después del final de rrnS en esa orientación.
# ------------------------------------------------------------

rrns_end_rc = L - start
cut_reverse = (rrns_end_rc + OFFSET_AFTER_RRNS) % L

reverse_ok, reverse_order = canonical_order_after_cut(
    cut_reverse,
    marker_features,
    L,
    orientation="reverse"
)

# ------------------------------------------------------------
# Elegir orientación que produzca:
# región control -> trnI -> trnQ -> trnM -> ...
# ------------------------------------------------------------

if forward_ok and not reverse_ok:
    rotated_seq = seq_to_use[cut_forward:] + seq_to_use[:cut_forward]
    cut_pos = cut_forward
    operation = "rotate_rrnS_canonical_forward"

elif reverse_ok and not forward_ok:
    rc_seq = seq_to_use.reverse_complement()
    rotated_seq = rc_seq[cut_reverse:] + rc_seq[:cut_reverse]
    cut_pos = cut_reverse
    operation = "reverse_complement_then_rotate_rrnS_canonical"

elif reverse_ok and forward_ok:
    # Caso raro: ambas orientaciones parecen compatibles.
    # Se elige la que deja trnI más cerca después del corte.
    if reverse_order["trnI"] < forward_order["trnI"]:
        rc_seq = seq_to_use.reverse_complement()
        rotated_seq = rc_seq[cut_reverse:] + rc_seq[:cut_reverse]
        cut_pos = cut_reverse
        operation = "reverse_complement_then_rotate_rrnS_canonical_both_possible"
    else:
        rotated_seq = seq_to_use[cut_forward:] + seq_to_use[:cut_forward]
        cut_pos = cut_forward
        operation = "rotate_rrnS_canonical_forward_both_possible"

else:
    # Fallback conservador:
    # Si no se pudieron validar trnI/trnQ/trnM, usar la hebra anotada de rrnS.
    # Esto evita romper muestras con anotaciones incompletas.
    if strand == -1:
        rc_seq = seq_to_use.reverse_complement()
        rotated_seq = rc_seq[cut_reverse:] + rc_seq[:cut_reverse]
        cut_pos = cut_reverse
        operation = "reverse_complement_then_rotate_rrnS_fallback_by_rrnS_strand"
    else:
        rotated_seq = seq_to_use[cut_forward:] + seq_to_use[:cut_forward]
        cut_pos = cut_forward
        operation = "rotate_rrnS_fallback_by_rrnS_strand"

out_rec = SeqRecord(
    rotated_seq,
    id=f"{sample}_rrnS",
    name=f"{sample}_rrnS",
    description=(
        f"oriented_to_start_10nt_after_rrnS_end "
        f"expected_order=control_region-trnI-trnQ-trnM "
        f"rrnS_start_1based={start+1} "
        f"rrnS_end_1based={end} "
        f"strand={strand} "
        f"cut_position_1based={cut_pos+1} "
        f"operation={operation} "
        f"source={sequence_used}"
    )
)

output_fasta.parent.mkdir(parents=True, exist_ok=True)
SeqIO.write(out_rec, str(output_fasta), "fasta")

print(
    f"{sample}\t{genus}\t{input_fasta}\t{gb_used}\t"
    f"{start+1}\t{end}\t{strand}\t{cut_pos+1}\t"
    f"{len(input_seq)}\t{len(gb_seq)}\t"
    f"{sequence_used}\t{operation}\t{output_fasta}\tROTATED_OK"
)
PY
        then
          echo "  ✓ FASTA orientado 10 nt después de rrnS con orden canónico"

          if tail -n 1 "$REPORT" | grep -q "reverse_complement_then_rotate"; then
            echo "  ✓ Barcode invertido corregido con reverse-complement"
            ((reverse_ok+=1))
          fi

          ((rotated_ok+=1))
        else
          echo "  ! No se pudo ubicar rrnS o validar orientación. Se copiará sin rotar."
          cp "$CONSENSUS" "$RRNS_FASTA"
          echo -e "${SAMPLE_NAME}\t${GENUS}\t${CONSENSUS}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tinput_fasta_unrotated\tcopy_unrotated\t${RRNS_FASTA}\tRRNS_NOT_FOUND_OR_FAILED" >> "$REPORT"
          ((rrns_not_found+=1))
        fi
      fi
    fi

    if [[ ! -s "$RRNS_FASTA" ]]; then
      echo "  ! No se generó FASTA válido"
      echo
      continue
    fi

    if [[ "$RUN_MITOFINDER" == "1" ]]; then
      MF_RRNS_OUTDIR="$OUT_MITOFINDER_RRNS/$GENUS/$SAMPLE_NAME"
      mkdir -p "$MF_RRNS_OUTDIR"

      SAMPLE_RRNS_ID="${SAMPLE_SAFE}_rrnS"

      echo "  → Ejecutando MitoFinder sobre FASTA orientado"
      echo "    Input:  $RRNS_FASTA"
      echo "    Output: $MF_RRNS_OUTDIR"

      if ! (
        cd "$MF_RRNS_OUTDIR"

        mitofinder \
          --seqid "$SAMPLE_RRNS_ID" \
          --assembly "$RRNS_FASTA" \
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

      COI_FASTA="$MF_RRNS_OUTDIR/COI_${SAMPLE_RRNS_ID}.fna"

      mapfile -t CDS_FILES < <(
        find "$MF_RRNS_OUTDIR" -type f -path "*/CDS/*.fna" 2>/dev/null | sort
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
          echo "  ! No se encontró COI/COX1 en la corrida rrnS"
        fi
      else
        echo "  ! No se encontraron CDS/*.fna en la corrida rrnS"
      fi
    fi

    ((processed+=1))
    echo

  done
done

echo "============================================================"
echo "PROCESO COMPLETADO"
echo "Muestras procesadas:          $processed"
echo "Reordenadas por rrnS:         $rotated_ok"
echo "Invertidas corregidas:        $reverse_ok"
echo "rrnS no encontrado/fallido:   $rrns_not_found"
echo "Consensus faltantes:          $missing_consensus"
echo "GenBank previos faltantes:    $missing_genbank"
echo "MitoFinder fallidos:          $failed_mitofinder"
echo
echo "FASTA orientados en:"
echo "  $OUT_RRNS"
echo
echo "MitoFinder rrnS en:"
echo "  $OUT_MITOFINDER_RRNS"
echo
echo "Reporte:"
echo "  $REPORT"
echo "============================================================"
