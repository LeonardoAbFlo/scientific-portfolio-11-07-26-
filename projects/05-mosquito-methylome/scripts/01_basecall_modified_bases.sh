#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

usage() {
  cat <<'EOF'
Usage: 01_basecall_modified_bases.sh --pod5 DIR --output DIR --model MODEL --kit KIT [options]

Options:
  --device DEVICE       Dorado device (default: auto)
  --min-qscore NUMBER   Minimum Q score (default: 10)
  --force               Replace existing BAM files
EOF
}

pod5_dir=""; output=""; model=""; kit=""; device="auto"; min_qscore=10; force=0
while (($#)); do
  case "$1" in
    --pod5) pod5_dir="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --model) model="$2"; shift 2 ;;
    --kit) kit="$2"; shift 2 ;;
    --device) device="$2"; shift 2 ;;
    --min-qscore) min_qscore="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done
[[ -n "$pod5_dir" && -n "$output" && -n "$model" && -n "$kit" ]] || { usage >&2; exit 2; }
require_directory "$pod5_dir"
require_command dorado
mkdir -p "$output"
init_log "$output/01_basecall.log"

pod5_files=("$pod5_dir"/*.pod5)
((${#pod5_files[@]})) || die "No POD5 files found in $pod5_dir"
for pod5 in "${pod5_files[@]}"; do
  output_bam="$output/$(basename "${pod5%.pod5}").bam"
  if [[ -s "$output_bam" && $force -eq 0 ]]; then
    log "SKIP $(basename "$pod5")"
    continue
  fi
  temp_bam="${output_bam}.partial"
  rm -f -- "$temp_bam"
  log "Basecalling $(basename "$pod5")"
  if dorado basecaller "$model" "$pod5" --kit-name "$kit" --device "$device" --min-qscore "$min_qscore" > "$temp_bam"; then
    mv "$temp_bam" "$output_bam"
  else
    rm -f -- "$temp_bam"
    die "Dorado failed for $pod5"
  fi
done
dorado --version > "$output/dorado_version.txt" 2>&1 || true
log "Basecalling complete"

