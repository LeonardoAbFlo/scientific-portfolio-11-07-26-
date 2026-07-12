#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

usage() { printf 'Usage: 05_call_methylated_sites.sh --alignments DIR --output DIR [--threads N] [--filter FRACTION] [--force]\n'; }
alignments=""; output=""; threads=16; filter=0.5; force=0
while (($#)); do
  case "$1" in
    --alignments) alignments="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --filter) filter="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done
[[ -n "$alignments" && -n "$output" ]] || { usage >&2; exit 2; }
require_directory "$alignments"; require_command modkit
mkdir -p "$output"
init_log "$output/05_call_methylated_sites.log"

bams=("$alignments"/*.bam)
((${#bams[@]})) || die "No BAM files found in $alignments"
for bam in "${bams[@]}"; do
  sample="$(basename "$bam" .bam)"
  bed="$output/$sample.bedmethyl"
  if [[ -s "$bed" && $force -eq 0 ]]; then log "SKIP $sample"; continue; fi
  log "Modkit pileup: $sample"
  modkit pileup --filter-threshold "$filter" --threads "$threads" "$bam" "$bed"
done
modkit --version > "$output/modkit_version.txt" 2>&1 || true
log "Methylation calling complete"

