#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

usage() { printf 'Usage: 04_align_reads_to_mitogenomes.sh --reads DIR --samples TSV --output DIR [--threads N] [--force]\n'; }
reads=""; samples=""; output=""; threads=16; force=0
while (($#)); do
  case "$1" in
    --reads) reads="$2"; shift 2 ;;
    --samples) samples="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done
[[ -n "$reads" && -n "$samples" && -n "$output" ]] || { usage >&2; exit 2; }
require_directory "$reads"; require_file "$samples"
require_command minimap2; require_command samtools
mkdir -p "$output"
init_log "$output/04_align_reads.log"

while IFS=$'\t' read -r sample_id species barcode reference; do
  require_file "$reference"
  mapfile -t matches < <(find "$reads" -maxdepth 1 -type f \( -name "*${barcode}*.fastq" -o -name "*${barcode}*.fastq.gz" \) -print | sort)
  ((${#matches[@]})) || die "No reads matching $barcode for $sample_id"
  bam="$output/$sample_id.bam"
  if [[ -s "$bam" && -s "$bam.bai" && $force -eq 0 ]]; then log "SKIP $sample_id"; continue; fi
  temp="${bam}.partial"
  log "Aligning $sample_id ($species)"
  minimap2 -t "$threads" -a -y -x map-ont "$reference" "${matches[@]}" \
    | samtools sort -@ "$threads" -o "$temp" -
  mv "$temp" "$bam"
  samtools index -@ "$threads" "$bam"
done < <(sample_rows "$samples")
log "Alignment complete"

