#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

usage() {
  cat <<'EOF'
Usage: 04_classify_contigs.sh --contigs DIR --output DIR [options]

Options:
  --threads N       Worker threads (default: 16)
  --tools LIST      Comma-separated subset of gtdbtk,checkm2 (default: both)
  --gtdb-data DIR   GTDB-Tk database; otherwise use GTDBTK_DATA_PATH
  --force           Replace existing tool output directories
EOF
}

contigs=""
output=""
threads=16
tools="gtdbtk,checkm2"
gtdb_data="${GTDBTK_DATA_PATH:-}"
force=0
while (($#)); do
  case "$1" in
    --contigs) contigs="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --tools) tools="$2"; shift 2 ;;
    --gtdb-data) gtdb_data="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$contigs" && -n "$output" ]] || { usage >&2; exit 2; }
require_directory "$contigs"
mkdir -p "$output"
init_log "$output/04_classify_contigs.log"

run_gtdb=0
run_checkm=0
IFS=',' read -r -a selected_tools <<< "$tools"
for tool in "${selected_tools[@]}"; do
  case "$tool" in
    gtdbtk) run_gtdb=1; require_command gtdbtk; [[ -n "$gtdb_data" ]] || die "Set --gtdb-data or GTDBTK_DATA_PATH" ;;
    checkm2) run_checkm=1; require_command checkm2 ;;
    *) die "Unsupported tool: $tool" ;;
  esac
done

export GTDBTK_DATA_PATH="$gtdb_data"
for sample_dir in "$contigs"/*; do
  [[ -d "$sample_dir" ]] || continue
  compgen -G "$sample_dir/*.fasta" >/dev/null || continue
  sample="$(basename "$sample_dir")"

  if ((run_gtdb)); then
    gtdb_out="$output/gtdbtk/$sample"
    if [[ ! -d "$gtdb_out" || $force -eq 1 ]]; then
      ((force)) && rm -rf -- "$gtdb_out"
      log "GTDB-Tk: $sample"
      gtdbtk classify_wf --genome_dir "$sample_dir" --out_dir "$gtdb_out" --cpus "$threads" -x fasta
    else
      log "SKIP GTDB-Tk: $sample"
    fi
  fi

  if ((run_checkm)); then
    checkm_out="$output/checkm2/$sample"
    if [[ ! -s "$checkm_out/quality_report.tsv" || $force -eq 1 ]]; then
      ((force)) && rm -rf -- "$checkm_out"
      log "CheckM2: $sample"
      checkm2 predict --input "$sample_dir" --output-directory "$checkm_out" --threads "$threads" -x fasta
    else
      log "SKIP CheckM2: $sample"
    fi
  fi
done

{
  ((run_gtdb)) && gtdbtk --version 2>&1 || true
  ((run_checkm)) && checkm2 --version 2>&1 || true
} > "$output/software_versions.txt"
log "Contig classification and quality assessment complete"

