#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

usage() { printf 'Usage: 02_demultiplex_reads.sh --basecalls DIR --output DIR [--force]\n'; }
basecalls=""; output=""; force=0
while (($#)); do
  case "$1" in
    --basecalls) basecalls="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done
[[ -n "$basecalls" && -n "$output" ]] || { usage >&2; exit 2; }
require_directory "$basecalls"
require_command dorado
if [[ -d "$output" && $force -eq 1 ]]; then rm -rf -- "$output"; fi
mkdir -p "$output"
init_log "$output/02_demultiplex.log"
log "Demultiplexing basecalls"
dorado demux --no-classify --output-dir "$output" "$basecalls"
log "Demultiplexing complete"

