#!/usr/bin/env bash
# Basecall POD5 files into one uBAM per input file.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$SCRIPT_DIR}"

FLOWCELL="${FLOWCELL:-F1}"
RUN_DIR="${RUN_DIR:-$WORKSPACE_ROOT/runs/$FLOWCELL/raw}"
OUT_DIR="${OUT_DIR:-$WORKSPACE_ROOT/runs/$FLOWCELL/basecalls}"
MODEL="sup@v5.0.0,6mA@v3"
KIT="SQK-NBD114-96"
MIN_Q=10
DEVICE="cuda:0,1"
LOG_DIR="${LOG_DIR:-$WORKSPACE_ROOT/logs}"
LOG_PATH="${LOG_DIR}/01_basecall_reads.${FLOWCELL}.log"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dorado

mkdir -p "$OUT_DIR" "$LOG_DIR"
: > "$LOG_PATH"

processed=0
skipped=0
failed=0

echo "[START] Basecalling flowcell ${FLOWCELL}" | tee -a "$LOG_PATH"

for pod5 in "$RUN_DIR"/*.pod5; do
    pod5_name="$(basename "$pod5")"
    bam_out="${OUT_DIR}/${pod5_name%.pod5}.bam"

    if [[ -e "$bam_out" ]]; then
        echo "[SKIP] ${pod5_name} already has a uBAM" | tee -a "$LOG_PATH"
        ((skipped+=1))
        continue
    fi

    echo "[INFO] Basecalling ${pod5_name}" | tee -a "$LOG_PATH"
    if dorado basecaller "$MODEL" "$pod5" \
           --kit-name "$KIT" \
           --device "$DEVICE" \
           --min-qscore "$MIN_Q" \
           > "$bam_out" 2>> "$LOG_PATH"; then
        echo "[OK] ${pod5_name}" | tee -a "$LOG_PATH"
        ((processed+=1))
    else
        echo "[FAIL] ${pod5_name}" | tee -a "$LOG_PATH"
        rm -f "$bam_out"
        ((failed+=1))
    fi
    echo "---" >> "$LOG_PATH"
done

echo
echo "========================================"
echo "[DONE] Basecalling complete"
echo "Flowcell: $FLOWCELL"
echo "Processed: $processed"
echo "Skipped: $skipped"
echo "Failed: $failed"
echo "Output: $OUT_DIR"
echo "Log: $LOG_PATH"
echo "========================================"
