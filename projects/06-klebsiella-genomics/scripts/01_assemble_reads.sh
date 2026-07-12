#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

usage() {
  cat <<'EOF'
Usage: 01_assemble_reads.sh --reads DIR --output DIR [options]

Options:
  --threads N       Flye threads (default: 16)
  --mode MODE       nano-hq, nano-raw, or nano-corr (default: nano-hq)
  --meta            Enable Flye metagenome mode
  --force           Replace an existing per-sample result
EOF
}

reads_dir=""
output_dir=""
threads=16
mode="nano-hq"
meta=0
force=0

while (($#)); do
  case "$1" in
    --reads) reads_dir="$2"; shift 2 ;;
    --output) output_dir="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --mode) mode="$2"; shift 2 ;;
    --meta) meta=1; shift ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$reads_dir" && -n "$output_dir" ]] || { usage >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
[[ "$mode" =~ ^nano-(hq|raw|corr)$ ]] || die "Unsupported Flye mode: $mode"
require_directory "$reads_dir"
require_command flye

mkdir -p "$output_dir"
init_log "$output_dir/01_assemble_reads.log"

inputs=("$reads_dir"/*.fastq "$reads_dir"/*.fastq.gz "$reads_dir"/*.fq "$reads_dir"/*.fq.gz)
((${#inputs[@]})) || die "No FASTQ files found in $reads_dir"

extra=()
((meta)) && extra+=(--meta)

for reads in "${inputs[@]}"; do
  sample="$(sample_name_from_fastq "$reads")"
  sample_out="$output_dir/$sample"
  assembly="$sample_out/assembly.fasta"

  if [[ -s "$assembly" && $force -eq 0 ]]; then
    log "SKIP $sample: $assembly already exists"
    continue
  fi
  if [[ -e "$sample_out" && $force -eq 1 ]]; then
    rm -rf -- "$sample_out"
  fi

  log "Assembling $sample"
  flye "--$mode" "$reads" --out-dir "$sample_out" --threads "$threads" "${extra[@]}"
done

flye --version > "$output_dir/software_versions.txt" 2>&1 || true
log "Assembly stage complete"

