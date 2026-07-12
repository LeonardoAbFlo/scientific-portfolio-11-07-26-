#!/usr/bin/env bash
# Run NanoPlot for each barcode FASTQ in the workspace.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$SCRIPT_DIR}"
RUNS_ROOT="${RUNS_ROOT:-$WORKSPACE_ROOT/runs}"
LOG_DIR="${LOG_DIR:-$WORKSPACE_ROOT/logs}"
LOG_PATH="${LOG_DIR}/04_assess_read_quality.log"

THREADS=30
CONDA_SH="/path/to/anaconda3/etc/profile.d/conda.sh"
CONDA_ENV="nanoplot_env"

source "$CONDA_SH"
conda activate "$CONDA_ENV"

mkdir -p "$LOG_DIR"
exec > >(tee -a "$LOG_PATH") 2>&1

echo "[INFO] Using NanoPlot: $(NanoPlot --version 2>/dev/null || true)"
echo "[INFO] Runs root: $RUNS_ROOT"
echo "[INFO] Threads: $THREADS"

processed_runs=0
processed_fastqs=0

for FC_DIR in "$RUNS_ROOT"/*; do
  [[ -d "$FC_DIR" ]] || continue
  FC="$(basename "$FC_DIR")"
  FASTQS=()

  while IFS= read -r -d '' f; do FASTQS+=("$f"); done < <(
    find "$FC_DIR/reads" -maxdepth 1 -type f -name "barcode*.fastq.gz" -print0 2>/dev/null || true
  )

  while IFS= read -r -d '' f; do FASTQS+=("$f"); done < <(
    find "$FC_DIR/reads" -mindepth 2 -maxdepth 2 -type f -name "barcode*.fastq.gz" -print0 2>/dev/null || true
  )

  if [[ ${#FASTQS[@]} -eq 0 ]]; then
    echo "[WARN] No FASTQ files found in $FC_DIR/reads"
    continue
  fi

  ((processed_runs+=1))
  echo "[INFO] Processing $FC with ${#FASTQS[@]} FASTQ files"

  for FQ in "${FASTQS[@]}"; do
    BC="$(basename "$FQ")"
    BC="${BC%.fastq.gz}"
    OUT="$FC_DIR/qc/$BC"
    mkdir -p "$OUT"

    echo "[INFO] Running NanoPlot for $FC/$BC"

    NanoPlot --fastq "$FQ" -o "$OUT" -t "$THREADS" --N50 --loglength \
      > "$OUT/nanoplot.log" 2>&1

    # Save a compact stats snapshot next to the NanoPlot output.
    {
      echo "sample=${FC}${BC}"
      echo "fastq=$FQ"
      echo "outdir=$OUT"
      echo "---- key stats ----"
      grep -i -E "N50|Mean read length|Median read length|Mean read quality|Median read quality|Reads" "$OUT"/* 2>/dev/null || true
    } > "$OUT/key_stats.txt"

    ((processed_fastqs+=1))
  done
done

echo
echo "========================================"
echo "[DONE] NanoPlot analysis complete"
echo "Runs: $processed_runs"
echo "FASTQ files: $processed_fastqs"
echo "Output: $RUNS_ROOT/*/qc"
echo "Log: $LOG_PATH"
echo "========================================"
