#!/usr/bin/env bash

die() {
  printf '[ERROR] %s\n' "$*" >&2
  exit 1
}

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

require_file() {
  [[ -f "$1" ]] || die "File not found: $1"
}

require_directory() {
  [[ -d "$1" ]] || die "Directory not found: $1"
}

init_log() {
  local log_file="$1"
  mkdir -p "$(dirname "$log_file")"
  exec > >(tee -a "$log_file") 2>&1
}

sample_rows() {
  local sample_sheet="$1"
  awk -F '\t' 'NR > 1 && $1 !~ /^#/ && NF >= 4 {print $1 "\t" $2 "\t" $3 "\t" $4}' "$sample_sheet"
}

