#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

usage() {
  cat <<'EOF'
Usage: 08_profile_filtered_assemblies.sh --assemblies DIR --output DIR [options]

Options:
  --threads N       Worker threads (default: 16)
  --tools LIST      Comma-separated mobsuite,rgi,checkm2,kleborate (default: all)
  --force           Replace existing results
EOF
}

assemblies=""
output=""
threads=16
tools="mobsuite,rgi,checkm2,kleborate"
force=0
while (($#)); do
  case "$1" in
    --assemblies) assemblies="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --tools) tools="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$assemblies" && -n "$output" ]] || { usage >&2; exit 2; }
require_directory "$assemblies"
mkdir -p "$output"
init_log "$output/08_profile_filtered_assemblies.log"

run_mob=0; run_rgi=0; run_checkm=0; run_kleborate=0
IFS=',' read -r -a selected_tools <<< "$tools"
for tool in "${selected_tools[@]}"; do
  case "$tool" in
    mobsuite) run_mob=1 ;;
    rgi) run_rgi=1 ;;
    checkm2) run_checkm=1; require_command checkm2 ;;
    kleborate) run_kleborate=1; require_command kleborate ;;
    *) die "Unsupported tool: $tool" ;;
  esac
done

profile_tools=()
((run_mob)) && profile_tools+=(mobsuite)
((run_rgi)) && profile_tools+=(rgi)
if ((${#profile_tools[@]})); then
  joined="$(IFS=,; printf '%s' "${profile_tools[*]}")"
  args=(--assemblies "$assemblies" --output "$output" --threads "$threads" --tools "$joined")
  ((force)) && args+=(--force)
  "$SCRIPT_DIR/03_profile_assemblies.sh" "${args[@]}"
fi

for sample_dir in "$assemblies"/*; do
  [[ -d "$sample_dir" && -s "$sample_dir/assembly.fasta" ]] || continue
  sample="$(basename "$sample_dir")"
  assembly="$sample_dir/assembly.fasta"

  if ((run_checkm)); then
    checkm_out="$output/checkm2/$sample"
    if [[ ! -s "$checkm_out/quality_report.tsv" || $force -eq 1 ]]; then
      ((force)) && rm -rf -- "$checkm_out"
      log "CheckM2: $sample"
      checkm2 predict --input "$sample_dir" --output-directory "$checkm_out" --threads "$threads" -x fasta
    fi
  fi

  if ((run_kleborate)); then
    kleborate_out="$output/kleborate/$sample.tsv"
    if [[ ! -s "$kleborate_out" || $force -eq 1 ]]; then
      mkdir -p "$(dirname "$kleborate_out")"
      log "Kleborate: $sample"
      kleborate -a "$assembly" -o "$kleborate_out" -p kpsc --trim_headers
    fi
  fi
done

log "Filtered assembly profiling complete"
