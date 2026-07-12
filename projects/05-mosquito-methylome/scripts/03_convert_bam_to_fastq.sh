#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

usage() { printf 'Usage: 03_convert_bam_to_fastq.sh --bams DIR --output DIR [--threads N] [--force]\n'; }
bams=""; output=""; threads=16; force=0
while (($#)); do
  case "$1" in
    --bams) bams="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done
[[ -n "$bams" && -n "$output" ]] || { usage >&2; exit 2; }
require_directory "$bams"
require_command samtools
mkdir -p "$output"
init_log "$output/03_convert_bam_to_fastq.log"

if command -v pigz >/dev/null 2>&1; then compressor=(pigz -p "$threads"); else compressor=(gzip); fi
bam_files=("$bams"/*.bam)
((${#bam_files[@]})) || die "No BAM files found in $bams"
for bam in "${bam_files[@]}"; do
  name="$(basename "$bam" .bam)"
  destination="$output/$name.fastq.gz"
  if [[ -s "$destination" && $force -eq 0 ]]; then log "SKIP $name"; continue; fi
  temp="${destination}.partial"
  log "Extracting FASTQ: $name"
  samtools fastq -@ "$threads" -n -T '*' "$bam" | "${compressor[@]}" > "$temp"
  mv "$temp" "$destination"
done
log "FASTQ extraction complete"

