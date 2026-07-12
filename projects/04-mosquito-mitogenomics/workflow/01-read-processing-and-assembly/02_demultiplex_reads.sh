#!/usr/bin/env bash
# Demultiplex basecalled BAMs for each run in the workspace.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$SCRIPT_DIR}"
RUNS_ROOT="${RUNS_ROOT:-$WORKSPACE_ROOT/runs}"
LOG_DIR="${LOG_DIR:-$WORKSPACE_ROOT/logs}"
LOG_PATH="${LOG_DIR}/02_demultiplex_reads.log"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dorado

mkdir -p "$LOG_DIR"
exec > >(tee -a "$LOG_PATH") 2>&1

processed_runs=0

demux_all_bams() {
    local in_dir="$1"
    local out_dir="$2"
    mkdir -p "$out_dir"

    echo "[INFO] Demultiplexing $in_dir"
    dorado demux --no-classify --output-dir "$out_dir" "$in_dir"
    echo "[OK] Wrote demultiplexed BAMs to $out_dir"
}

for run_dir in "$RUNS_ROOT"/*; do
    [[ -d "$run_dir" ]] || continue

    in_dir="${run_dir}/basecalls"
    out_dir="${run_dir}/demultiplexed"
    [[ -d "$in_dir" ]] || continue

    demux_all_bams "$in_dir" "$out_dir"
    ((processed_runs+=1))
done

echo
echo "========================================"
echo "[DONE] Demultiplexing complete"
echo "Runs: $processed_runs"
echo "Output: $RUNS_ROOT/*/demultiplexed"
echo "Log: $LOG_PATH"
echo "========================================"
