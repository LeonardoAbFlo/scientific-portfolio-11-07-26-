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

require_directory() {
  [[ -d "$1" ]] || die "Directory not found: $1"
}

init_log() {
  local log_file="$1"
  mkdir -p "$(dirname "$log_file")"
  exec > >(tee -a "$log_file") 2>&1
}

sample_name_from_fastq() {
  local name
  name="$(basename "$1")"
  name="${name%.gz}"
  name="${name%.fastq}"
  name="${name%.fq}"
  printf '%s\n' "$name"
}

