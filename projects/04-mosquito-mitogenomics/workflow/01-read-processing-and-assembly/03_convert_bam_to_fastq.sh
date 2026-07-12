#!/usr/bin/env bash
# Convert demultiplexed BAM files into compressed FASTQ files.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$SCRIPT_DIR}"
RUNS_ROOT="${RUNS_ROOT:-$WORKSPACE_ROOT/runs}"
LOG_DIR="${LOG_DIR:-$WORKSPACE_ROOT/logs}"
LOG_PATH="${LOG_DIR}/03_convert_bam_to_fastq.log"

CONDA_SH="/path/to/anaconda3/etc/profile.d/conda.sh"
ENV="samtools"
THREADS=30

source "$CONDA_SH"
conda activate "$ENV"
mkdir -p "$LOG_DIR"
exec > >(tee -a "$LOG_PATH") 2>&1
shopt -s nullglob

if command -v pigz >/dev/null 2>&1; then
  COMPRESS=(pigz -p "$THREADS")
else
  COMPRESS=(gzip)
fi

processed_runs=0
written_fastqs=0

for RUN_DIR in "$RUNS_ROOT"/*; do
  [[ -d "$RUN_DIR" ]] || continue

  IN_DIR="${RUN_DIR}/demultiplexed"
  OUT_DIR="${RUN_DIR}/reads"
  [[ -d "$IN_DIR" ]] || continue

  mkdir -p "$OUT_DIR"
  n_found=0

  for BAM in "$IN_DIR"/*.bam; do
    ((n_found+=1))
    base="$(basename "$BAM" .bam)"
    OUT="$OUT_DIR/${base}.fastq.gz"

    echo "[INFO] Converting $BAM -> $OUT"
    # Keep all AUX tags in the FASTQ header comment.
    samtools fastq -@ "$THREADS" -n -T '*' "$BAM" | "${COMPRESS[@]}" > "$OUT"
    ((written_fastqs+=1))
  done

  if (( n_found == 0 )); then
    echo "[WARN] No BAM files found in $IN_DIR"
  else
    ((processed_runs+=1))
  fi
done

echo
echo "========================================"
echo "[DONE] FASTQ conversion complete"
echo "Runs: $processed_runs"
echo "FASTQ files: $written_fastqs"
echo "Output: $RUNS_ROOT/*/reads"
echo "Log: $LOG_PATH"
echo "========================================"
