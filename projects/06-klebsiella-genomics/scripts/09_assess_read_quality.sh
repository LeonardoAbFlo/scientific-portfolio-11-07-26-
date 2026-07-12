#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

usage() {
  printf 'Usage: 09_assess_read_quality.sh --reads DIR --output DIR [--threads N]\n'
}

reads=""; output=""; threads=16
while (($#)); do
  case "$1" in
    --reads) reads="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done
[[ -n "$reads" && -n "$output" ]] || { usage >&2; exit 2; }
require_directory "$reads"
require_command NanoStat
mkdir -p "$output/stats"
init_log "$output/09_assess_read_quality.log"

printf 'sample_id\tmean_read_length\tmean_read_quality\tmedian_read_length\tmedian_read_quality\tn50\n' > "$output/read_qc_summary.tsv"
inputs=("$reads"/*.fastq "$reads"/*.fastq.gz "$reads"/*.fq "$reads"/*.fq.gz)
((${#inputs[@]})) || die "No FASTQ files found in $reads"

for fastq in "${inputs[@]}"; do
  sample="$(sample_name_from_fastq "$fastq")"
  stats="$output/stats/$sample.txt"
  log "NanoStat: $sample"
  NanoStat --fastq "$fastq" -t "$threads" > "$stats"
  value() { awk -F: -v key="$1" '$1 ~ key {sub(/^[[:space:]]+/, "", $2); gsub(/,/, "", $2); print $2; exit}' "$stats"; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$(value '^Mean read length')" "$(value '^Mean read quality')" \
    "$(value '^Median read length')" "$(value '^Median read quality')" "$(value '^Read length N50')" \
    >> "$output/read_qc_summary.tsv"
done
log "Read QC complete"
