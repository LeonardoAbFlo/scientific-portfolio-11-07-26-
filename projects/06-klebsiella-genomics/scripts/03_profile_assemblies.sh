#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

usage() {
  cat <<'EOF'
Usage: 03_profile_assemblies.sh --assemblies DIR --output DIR [options]

Options:
  --threads N          Worker threads (default: 16)
  --tools LIST         Comma-separated subset of mobsuite,rgi (default: both)
  --force              Allow tools to replace existing results
EOF
}

assemblies=""
output=""
threads=16
tools="mobsuite,rgi"
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
init_log "$output/03_profile_assemblies.log"

run_mob=0
run_rgi=0
IFS=',' read -r -a selected_tools <<< "$tools"
for tool in "${selected_tools[@]}"; do
  case "$tool" in
    mobsuite) run_mob=1; require_command mob_recon ;;
    rgi) run_rgi=1; require_command rgi ;;
    *) die "Unsupported tool: $tool" ;;
  esac
done

for sample_dir in "$assemblies"/*; do
  [[ -d "$sample_dir" ]] || continue
  sample="$(basename "$sample_dir")"
  assembly="$sample_dir/assembly.fasta"
  [[ -s "$assembly" ]] || assembly="$sample_dir/consensus.fasta"
  [[ -s "$assembly" ]] || { log "WARN no assembly for $sample"; continue; }

  if ((run_mob)); then
    mob_out="$output/mobsuite/$sample"
    if [[ ! -s "$mob_out/contig_report.txt" || $force -eq 1 ]]; then
      mkdir -p "$mob_out"
      log "MOB-suite: $sample"
      mob_args=(--infile "$assembly" --outdir "$mob_out" -n "$threads" -c)
      ((force)) && mob_args+=(--force)
      mob_recon "${mob_args[@]}"
    else
      log "SKIP MOB-suite: $sample"
    fi
  fi

  if ((run_rgi)); then
    rgi_out="$output/rgi/$sample"
    rgi_prefix="$rgi_out/$sample"
    if [[ ! -s "${rgi_prefix}.txt" || $force -eq 1 ]]; then
      mkdir -p "$rgi_out"
      log "RGI: $sample"
      rgi_args=(main --input_sequence "$assembly" --output_file "$rgi_prefix" --local -n "$threads")
      ((force)) && rgi_args+=(--clean)
      rgi "${rgi_args[@]}"
    else
      log "SKIP RGI: $sample"
    fi
  fi
done

{
  ((run_mob)) && mob_recon --version 2>&1 || true
  ((run_rgi)) && rgi --version 2>&1 || true
} > "$output/software_versions.txt"
log "Initial assembly profiling complete"

