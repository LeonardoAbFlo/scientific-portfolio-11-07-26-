#!/usr/bin/env bash
# Filter barcode FASTQs against barcode-specific mitochondrial references.
set -euo pipefail

THREADS=30

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$SCRIPT_DIR}"
RUNS_ROOT="${RUNS_ROOT:-$WORKSPACE_ROOT/runs}"

CONDA_ENV="minimap2"
CONDA_SH="/path/to/anaconda3/etc/profile.d/conda.sh"

declare -A FLOWCELL_BARCODES
FLOWCELL_BARCODES["F1"]="81 83 84"
FLOWCELL_BARCODES["F3"]="89 90 92 93 95"
FLOWCELL_BARCODES["F4"]="90 93"

LOG_DIR="${WORKSPACE_ROOT}/logs"
mkdir -p "$LOG_DIR"

SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
LOG_PATH="${LOG_DIR}/${SCRIPT_NAME%.sh}.log"

exec > >(tee -a "$LOG_PATH") 2>&1

echo "[START] Mitochondrial read filtering"
echo "[INFO] Workspace root: $WORKSPACE_ROOT"
echo "[INFO] Threads: $THREADS"
echo "[INFO] Conda env: $CONDA_ENV"

if [[ -f "$CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$CONDA_SH"
else
  echo "ERROR: Missing $CONDA_SH"
  exit 1
fi

conda activate "$CONDA_ENV" || {
  echo "ERROR: Cannot activate conda env $CONDA_ENV"
  exit 1
}

command -v minimap2 >/dev/null 2>&1 || {
  echo "ERROR: minimap2 is not available in $CONDA_ENV"
  exit 1
}

command -v samtools >/dev/null 2>&1 || {
  echo "ERROR: samtools is not available in $CONDA_ENV"
  exit 1
}

command -v gzip >/dev/null 2>&1 || {
  echo "ERROR: gzip is not available"
  exit 1
}

command -v zcat >/dev/null 2>&1 || {
  echo "ERROR: zcat is not available"
  exit 1
}

PROCESSED=()
ALREADY_FILTERED=()
SKIPPED_REFS=()
SKIPPED_FASTQS=()
EMPTY_OUTPUTS=()

GLOBAL_START=$(date +%s)

for FLOWCELL in "${!FLOWCELL_BARCODES[@]}"; do
  echo
  echo "[INFO] Flowcell: $FLOWCELL"

  FASTQ_DIR="${RUNS_ROOT}/${FLOWCELL}/reads"
  REF_DIR="${RUNS_ROOT}/${FLOWCELL}/references"
  OUT_DIR="${RUNS_ROOT}/${FLOWCELL}/filtered_reads"
  mkdir -p "$OUT_DIR"

  for BARCODE_NUM in ${FLOWCELL_BARCODES[$FLOWCELL]}; do
    BARCODE_NAME="barcode${BARCODE_NUM}"
    FASTQ_FILE="${FASTQ_DIR}/${BARCODE_NAME}.fastq.gz"
    OUT_FASTQ="${OUT_DIR}/${BARCODE_NAME}.fastq.gz"

    echo "[INFO] Sample: ${FLOWCELL}/${BARCODE_NAME}"
    echo "[INFO] Input FASTQ: $FASTQ_FILE"
    echo "[INFO] Output FASTQ: $OUT_FASTQ"

    if [[ -s "$OUT_FASTQ" ]]; then
      echo "[SKIP] Filtered FASTQ already exists"

      if gzip -t "$OUT_FASTQ" 2>/dev/null; then
        READ_COUNT=$(zcat "$OUT_FASTQ" | awk 'END {print NR/4}')
        echo "[INFO] Existing reads: $READ_COUNT"
        ALREADY_FILTERED+=("${FLOWCELL}/${BARCODE_NAME} reads=${READ_COUNT}")
      else
        echo "[WARN] Existing file is not a valid gzip archive: $OUT_FASTQ"
        ALREADY_FILTERED+=("${FLOWCELL}/${BARCODE_NAME} invalid_gzip")
      fi

      continue
    fi

    if [[ ! -f "$FASTQ_FILE" ]]; then
      echo "[WARN] FASTQ not found"
      SKIPPED_FASTQS+=("${FLOWCELL}/${BARCODE_NAME}")
      continue
    fi

    shopt -s nullglob
    REF_MATCHES=( "${REF_DIR}/b${BARCODE_NUM}_"*.fasta "${REF_DIR}/b${BARCODE_NUM}_"*.fa )
    shopt -u nullglob

    if [[ ${#REF_MATCHES[@]} -eq 0 ]]; then
      echo "[WARN] Reference not found"
      echo "[INFO] Expected pattern: ${REF_DIR}/b${BARCODE_NUM}_*.fasta"
      SKIPPED_REFS+=("${FLOWCELL}/${BARCODE_NAME}")
      continue
    fi

    if [[ ${#REF_MATCHES[@]} -gt 1 ]]; then
      echo "[WARN] Multiple references found, using the first match"
      printf '  %s\n' "${REF_MATCHES[@]}"
    fi

    REF_FASTA="${REF_MATCHES[0]}"
    echo "[INFO] Reference: $REF_FASTA"

    START=$(date +%s)

    minimap2 \
      -t "$THREADS" \
      -ax map-ont \
      "$REF_FASTA" \
      "$FASTQ_FILE" \
      | samtools view \
          -@ "$THREADS" \
          -b \
          -F 2308 \
          - \
      | samtools fastq \
          -@ "$THREADS" \
          -n \
          - \
      | gzip -c > "$OUT_FASTQ"

    END=$(date +%s)
    DURATION_MIN=$(awk "BEGIN {printf \"%.2f\", ($END - $START)/60}")

    if [[ -s "$OUT_FASTQ" ]]; then
      if gzip -t "$OUT_FASTQ" 2>/dev/null; then
        READ_COUNT=$(zcat "$OUT_FASTQ" | awk 'END {print NR/4}')
        echo "[OK] Finished in ${DURATION_MIN} minutes"
        echo "[INFO] Filtered reads: $READ_COUNT"
        PROCESSED+=("${FLOWCELL}/${BARCODE_NAME} reads=${READ_COUNT}")
      else
        echo "[WARN] Output FASTQ is not a valid gzip archive: $OUT_FASTQ"
        EMPTY_OUTPUTS+=("${FLOWCELL}/${BARCODE_NAME} invalid_gzip")
      fi
    else
      echo "[WARN] Output FASTQ is empty"
      EMPTY_OUTPUTS+=("${FLOWCELL}/${BARCODE_NAME} reads=0")
    fi
  done
done

GLOBAL_END=$(date +%s)
TOTAL_MIN=$(awk "BEGIN {printf \"%.2f\", ($GLOBAL_END - $GLOBAL_START)/60}")

echo
echo "========================================"
echo "[DONE] Mitochondrial filtering complete"
echo "Processed: ${#PROCESSED[@]}"
echo "Skipped existing: ${#ALREADY_FILTERED[@]}"
echo "Skipped missing FASTQ: ${#SKIPPED_FASTQS[@]}"
echo "Skipped missing reference: ${#SKIPPED_REFS[@]}"
echo "Invalid or empty outputs: ${#EMPTY_OUTPUTS[@]}"
echo "Runtime (min): $TOTAL_MIN"
echo "Output: ${RUNS_ROOT}/*/filtered_reads"
echo "Log: $LOG_PATH"
echo "========================================"
