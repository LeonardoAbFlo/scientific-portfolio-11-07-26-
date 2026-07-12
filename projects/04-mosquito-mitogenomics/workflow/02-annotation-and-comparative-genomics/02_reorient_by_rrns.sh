#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Rotate Medaka consensuses to the rrnS reference point using the first MitoFinder pass.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd -P)}"
ANNOTATION_ROOT="${ANNOTATION_ROOT:-$WORKSPACE_ROOT/results/annotation}"
REPORT_ROOT="${REPORT_ROOT:-$WORKSPACE_ROOT/results/reports/module2}"
LOG_DIR="${LOG_DIR:-$WORKSPACE_ROOT/logs}"
LOG_PATH="${LOG_DIR}/02_reorient_by_rrns.log"

INPUT_MEDAKA_ROOT="${INPUT_MEDAKA_ROOT:-$ANNOTATION_ROOT/mitogenomes_medaka}"
INPUT_MITOFINDER_ROOT="${INPUT_MITOFINDER_ROOT:-$ANNOTATION_ROOT/medaka_mitofinder}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$ANNOTATION_ROOT/mitogenomes_medaka_rrnS}"
REPORT="${REPORT:-$REPORT_ROOT/02_reorient_by_rrns.report.tsv}"

RRNS_OFFSET="${RRNS_OFFSET:-10}"
LINE_WIDTH="${LINE_WIDTH:-80}"
FORCE="${FORCE:-0}"

usage() {
  cat <<'EOF'
Usage:
  02_reorient_by_rrns.sh

Key variables:
  INPUT_MEDAKA_ROOT=/path/to/mitogenomes_medaka
  INPUT_MITOFINDER_ROOT=/path/to/medaka_mitofinder
  OUTPUT_ROOT=/path/to/mitogenomes_medaka_rrnS
  REPORT=/path/to/02_reorient_by_rrns.report.tsv
  RRNS_OFFSET=10
  LINE_WIDTH=80
  FORCE=1
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

mkdir -p "$OUTPUT_ROOT" "$(dirname "$REPORT")" "$LOG_DIR"
exec > >(tee -a "$LOG_PATH") 2>&1

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
  log "ERROR: $*"
  exit 1
}

read_fasta_sequence() {
  local fasta_path="$1"
  awk 'BEGIN { ORS="" } /^>/ { next } { gsub(/[[:space:]]/, "", $0); printf "%s", toupper($0) }' "$fasta_path"
}

write_fasta() {
  local sequence_id="$1"
  local sequence="$2"
  local output_path="$3"
  local index=0
  local seq_len="${#sequence}"

  {
    printf '>%s\n' "$sequence_id"
    while (( index < seq_len )); do
      printf '%s\n' "${sequence:index:LINE_WIDTH}"
      ((index += LINE_WIDTH))
    done
  } > "$output_path"
}

extract_rrns_feature() {
  local gff_path="$1"

  awk '
    BEGIN { IGNORECASE=1; FS=OFS="\t"; best_priority=999 }
    /^#/ { next }
    (($3 == "rRNA") || ($3 == "gene")) && $9 ~ /(rrnS|12S)/ {
      priority = ($3 == "rRNA") ? 1 : 2
      if (priority < best_priority) {
        best = $4 OFS $5 OFS $7 OFS $3
        best_priority = priority
      }
    }
    END {
      if (best_priority != 999) {
        print best
      }
    }
  ' "$gff_path"
}

is_current() {
  local output="$1"
  shift
  local input

  [[ -s "$output" ]] || return 1
  for input in "$@"; do
    [[ "$output" -nt "$input" ]] || return 1
  done
}

if [[ ! -d "$INPUT_MEDAKA_ROOT" ]]; then
  die "Input Medaka directory does not exist: $INPUT_MEDAKA_ROOT"
fi

if [[ ! -d "$INPUT_MITOFINDER_ROOT" ]]; then
  die "Input MitoFinder directory does not exist: $INPUT_MITOFINDER_ROOT"
fi

printf 'genus\tsample\tinput_consensus\tannotation_gff\trrns_start\trrns_end\trrns_strand\tsequence_length_bp\tcut_position_0based\toutput_fasta\tstatus\n' > "$REPORT"

log "Workspace root: $WORKSPACE_ROOT"
log "Input Medaka root: $INPUT_MEDAKA_ROOT"
log "Input MitoFinder root: $INPUT_MITOFINDER_ROOT"
log "Output root: $OUTPUT_ROOT"
log "rrnS offset: $RRNS_OFFSET"
log "Report: $REPORT"

processed=0
rotated=0
unchanged=0
skipped_current=0
missing_consensus=0
missing_annotation=0
missing_rrns=0

for GENUS_DIR in "$INPUT_MEDAKA_ROOT"/*/; do
  [[ -d "$GENUS_DIR" ]] || continue
  GENUS="$(basename "$GENUS_DIR")"

  for SAMPLE_DIR in "$GENUS_DIR"*/; do
    [[ -d "$SAMPLE_DIR" ]] || continue
    SAMPLE_NAME="$(basename "$SAMPLE_DIR")"

    CONSENSUS="$SAMPLE_DIR/consensus.fasta"
    if [[ ! -s "$CONSENSUS" ]]; then
      log "WARN: Missing consensus FASTA for $GENUS/$SAMPLE_NAME"
      printf '%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_consensus\n' \
        "$GENUS" "$SAMPLE_NAME" "$CONSENSUS" >> "$REPORT"
      ((missing_consensus += 1))
      continue
    fi

    MF_SAMPLE_DIR="$INPUT_MITOFINDER_ROOT/$GENUS/$SAMPLE_NAME"
    if [[ ! -d "$MF_SAMPLE_DIR" ]]; then
      log "WARN: Missing MitoFinder directory for $GENUS/$SAMPLE_NAME"
      printf '%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_annotation\n' \
        "$GENUS" "$SAMPLE_NAME" "$CONSENSUS" >> "$REPORT"
      ((missing_annotation += 1))
      continue
    fi

    mapfile -t GFF_FILES < <(find "$MF_SAMPLE_DIR" -type f -name '*_mtDNA_contig.gff' | sort)

    if ((${#GFF_FILES[@]} == 0)); then
      log "WARN: Missing MitoFinder GFF for $GENUS/$SAMPLE_NAME"
      printf '%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_annotation\n' \
        "$GENUS" "$SAMPLE_NAME" "$CONSENSUS" >> "$REPORT"
      ((missing_annotation += 1))
      continue
    fi

    if ((${#GFF_FILES[@]} > 1)); then
      log "WARN: Multiple GFF files found for $GENUS/$SAMPLE_NAME; using the first one"
      printf '  %s\n' "${GFF_FILES[@]}"
    fi

    GFF_PATH="${GFF_FILES[0]}"
    RRNS_INFO="$(extract_rrns_feature "$GFF_PATH")"

    if [[ -z "$RRNS_INFO" ]]; then
      log "WARN: rrnS was not found in $GFF_PATH"
      printf '%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_rrnS\n' \
        "$GENUS" "$SAMPLE_NAME" "$CONSENSUS" "$GFF_PATH" >> "$REPORT"
      ((missing_rrns += 1))
      continue
    fi

    IFS=$'\t' read -r RRNS_START RRNS_END RRNS_STRAND RRNS_FEATURE_TYPE <<< "$RRNS_INFO"

    SEQUENCE="$(read_fasta_sequence "$CONSENSUS")"
    SEQ_LEN="${#SEQUENCE}"
    [[ "$SEQ_LEN" -gt 0 ]] || die "Consensus FASTA is empty after parsing: $CONSENSUS"

    CUT_POSITION=$(( (RRNS_END + RRNS_OFFSET) % SEQ_LEN ))

    if [[ "$SAMPLE_NAME" =~ ^(F[0-9]+)_([0-9]+)- ]]; then
      OUTPUT_ID="${BASH_REMATCH[1]}_barcode${BASH_REMATCH[2]}_rrnS"
    else
      SAMPLE_SAFE="$(echo "$SAMPLE_NAME" | sed -E 's/[^A-Za-z0-9_]+/_/g; s/_+$//')"
      OUTPUT_ID="${SAMPLE_SAFE}_rrnS"
    fi

    OUT_DIR="$OUTPUT_ROOT/$GENUS/$SAMPLE_NAME"
    OUT_FASTA="$OUT_DIR/${OUTPUT_ID}.fasta"
    METADATA_PATH="$OUT_DIR/reorientation_metadata.txt"
    mkdir -p "$OUT_DIR"

    if [[ "$FORCE" != "1" ]] && is_current "$OUT_FASTA" "$CONSENSUS" "$GFF_PATH"; then
      log "Skipping current output for $GENUS/$SAMPLE_NAME"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tskipped_current\n' \
        "$GENUS" "$SAMPLE_NAME" "$CONSENSUS" "$GFF_PATH" "$RRNS_START" "$RRNS_END" "$RRNS_STRAND" "$SEQ_LEN" "$CUT_POSITION" "$OUT_FASTA" >> "$REPORT"
      ((skipped_current += 1))
      continue
    fi

    if (( CUT_POSITION == 0 )); then
      ROTATED_SEQUENCE="$SEQUENCE"
      STATUS="already_oriented"
      ((unchanged += 1))
    else
      ROTATED_SEQUENCE="${SEQUENCE:CUT_POSITION}${SEQUENCE:0:CUT_POSITION}"
      STATUS="rotated"
      ((rotated += 1))
    fi

    write_fasta "$OUTPUT_ID" "$ROTATED_SEQUENCE" "$OUT_FASTA"

    {
      echo "sample=${GENUS}/${SAMPLE_NAME}"
      echo "input_consensus=$CONSENSUS"
      echo "annotation_gff=$GFF_PATH"
      echo "rrns_feature_type=$RRNS_FEATURE_TYPE"
      echo "rrns_start=$RRNS_START"
      echo "rrns_end=$RRNS_END"
      echo "rrns_strand=$RRNS_STRAND"
      echo "sequence_length_bp=$SEQ_LEN"
      echo "rrns_offset_bp=$RRNS_OFFSET"
      echo "cut_position_0based=$CUT_POSITION"
      echo "output_fasta=$OUT_FASTA"
      echo "status=$STATUS"
      echo "log=$LOG_PATH"
    } > "$METADATA_PATH"

    if [[ "$RRNS_STRAND" == "+" ]]; then
      log "WARN: rrnS is annotated on the plus strand for $GENUS/$SAMPLE_NAME; sequence was rotated without reverse complementing."
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$GENUS" "$SAMPLE_NAME" "$CONSENSUS" "$GFF_PATH" "$RRNS_START" "$RRNS_END" "$RRNS_STRAND" "$SEQ_LEN" "$CUT_POSITION" "$OUT_FASTA" "$STATUS" >> "$REPORT"

    ((processed += 1))
    log "Ready: $GENUS/$SAMPLE_NAME -> $OUT_FASTA"
  done
done

echo
echo "========================================"
echo "[DONE] rrnS reorientation complete"
echo "Processed: $processed"
echo "Rotated: $rotated"
echo "Already oriented: $unchanged"
echo "Skipped current: $skipped_current"
echo "Missing consensus: $missing_consensus"
echo "Missing annotation: $missing_annotation"
echo "Missing rrnS: $missing_rrns"
echo "Output: $OUTPUT_ROOT"
echo "Report: $REPORT"
echo "Log: $LOG_PATH"
echo "========================================"
