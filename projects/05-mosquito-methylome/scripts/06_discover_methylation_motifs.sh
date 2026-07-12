#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

usage() { printf 'Usage: 06_discover_methylation_motifs.sh --pileups DIR --samples TSV --output DIR [--threads N] [--force]\n'; }
pileups=""; samples=""; output=""; threads=16; force=0
while (($#)); do
  case "$1" in
    --pileups) pileups="$2"; shift 2 ;;
    --samples) samples="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done
[[ -n "$pileups" && -n "$samples" && -n "$output" ]] || { usage >&2; exit 2; }
require_directory "$pileups"; require_file "$samples"; require_command modkit
mkdir -p "$output"
init_log "$output/06_discover_methylation_motifs.log"

while IFS=$'\t' read -r sample_id species barcode reference; do
  pileup="$pileups/$sample_id.bedmethyl"
  motif="$output/$sample_id.tsv"
  require_file "$pileup"; require_file "$reference"
  if [[ -s "$motif" && $force -eq 0 ]]; then log "SKIP $sample_id"; continue; fi
  log "Motif search: $sample_id ($species; $barcode)"
  modkit motif search --in-bedmethyl "$pileup" --ref "$reference" -o "$motif" --threads "$threads"
done < <(sample_rows "$samples")
log "Motif discovery complete"
