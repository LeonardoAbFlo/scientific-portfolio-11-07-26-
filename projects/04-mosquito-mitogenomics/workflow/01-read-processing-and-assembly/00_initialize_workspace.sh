#!/usr/bin/env bash
set -euo pipefail

DEFAULT_FLOWCELLS=(F1 F2 F3 F4 F5)
RUNS_ROOT="runs"
FLOWCELL_STAGE_DIRS=(
  raw
  basecalls
  demultiplexed
  reads
  qc
  references
  references/reoriented_mitogenomes
  filtered_reads
  assemblies
  alignments
  alignments/dnaapler
  polishing
  annotation
  annotation/modkit_dnaapler
  annotation/modkit_motifs
  reports
  reports/methylation
)
GLOBAL_DIR_SUFFIXES=(
  config
  downloads
  downloads/zenodo_fastq
  logs
  references
  resources/models
  resources/mitos2_refdata
  resources/phylogeny
  resources/phylogeny/reference_fastas
  resources/phylogeny/outgroup_fastas
  resources/phylogeny/reference_maps
  work/assemblies
  work/mapping
  work/polishing
  work/annotation
  results/qc
  results/assemblies/selected
  results/mapping
  results/polishing/dorado_medaka
  results/polishing/map60_dorado_medaka
  results/annotation
  results/annotation/data
  results/annotation/ref_gb
  results/annotation/mitogenomes_medaka
  results/annotation/medaka_mitofinder
  results/annotation/mitogenomes_medaka_rrnS
  results/annotation/medaka_mitos2_rrnS
  results/annotation/corrections/mitogenomes_medaka_rrnS
  results/annotation/corrections/medaka_mitofinder_rrnS
  results/annotation/plots/circos
  results/annotation/plots/synteny
  results/annotation/stats
  results/phylogeny
  results/phylogeny/reference_mitofinder
  results/phylogeny/mitofinder_results
  results/phylogeny/concat
  results/phylogeny/tmp
  results/phylogeny/tree
  results/phylogeny/itol
  results/reports
  results/reports/module2
  results/reports/module3
  results/reports/module4
  tmp
)
LEGACY_FLOWCELL_ALIASES=(
  "bam:basecalls"
  "demux_bam:demultiplexed"
  "fastq_concat:reads"
  "nanoplot_results:qc"
  "ref_2map:references"
  "mitogenomas_reorientados:references/reoriented_mitogenomes"
  "p1_filtered_fastq:filtered_reads"
  "aln_dnaapler:alignments/dnaapler"
  "modkit_dnaapler:annotation/modkit_dnaapler"
  "modkit_motifs:annotation/modkit_motifs"
)
LEGACY_ROOT_ALIASES=(
  "dorado_models:resources/models"
  "mitos2_refdata:resources/mitos2_refdata"
  "paper1_flye_coverage_grid_select:work/assemblies"
  "paper1_final_selected_flye_sets:results/assemblies/selected"
  "dorado_medaka_polish_output:results/polishing/dorado_medaka"
  "map60_dorado_medaka_polish_output:results/polishing/map60_dorado_medaka"
  "maptest:results/mapping"
  "mitogenomas_medaka_rrnS:results/annotation/mitogenomes_medaka_rrnS"
  "medaka_MITOS2_rrnS:results/annotation/medaka_mitos2_rrnS"
  "mitofinder_results_Jimena:results/phylogeny/mitofinder_results"
  "concat:results/phylogeny/concat"
  "Jimena:results/phylogeny/itol"
  "results_p1/data:results/annotation/data"
  "results_p1/ref_gb:results/annotation/ref_gb"
  "results_p1/mitogenomes_medaka:results/annotation/mitogenomes_medaka"
  "results_p1/medaka_mitofinder:results/annotation/medaka_mitofinder"
  "results_p1/CORRECTIONS_mitogenomas_medaka_rrnS:results/annotation/corrections/mitogenomes_medaka_rrnS"
  "results_p1/CORRECTIONS_mitogenomes_medaka_rrnS:results/annotation/corrections/mitogenomes_medaka_rrnS"
  "results_p1/CORRECTIONS_medaka_mitofinder_rrnS:results/annotation/corrections/medaka_mitofinder_rrnS"
)

ROOT_DIR=""
declare -a FLOWCELLS=()

usage() {
  cat <<'EOF'
Usage:
  workspace.sh --root /path/to/workspace [--flowcells F1,F2,F3]
  workspace.sh /path/to/workspace [F1 F2 F3]

Description:
  Create a standard sequencing-analysis workspace layout with legacy
  compatibility aliases for the current pipeline scripts.
  Existing directories are preserved.

Default flowcells:
  F1 F2 F3 F4 F5

Examples:
  ./00_initialize_workspace.sh --root /home/user/mosquitos
  ./00_initialize_workspace.sh /home/user/mosquitos F1 F3 F4
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

append_flowcells() {
  local raw_list="$1"
  local entry

  for entry in ${raw_list//,/ }; do
    [[ -n "$entry" ]] || continue
    FLOWCELLS+=("$entry")
  done
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --root)
        [[ $# -ge 2 ]] || die "Missing value for --root"
        ROOT_DIR="$2"
        shift 2
        ;;
      --flowcells)
        [[ $# -ge 2 ]] || die "Missing value for --flowcells"
        append_flowcells "$2"
        shift 2
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      -*)
        die "Unknown option: $1"
        ;;
      *)
        if [[ -z "$ROOT_DIR" ]]; then
          ROOT_DIR="$1"
        else
          FLOWCELLS+=("$1")
        fi
        shift
        ;;
    esac
  done
}

normalize_flowcells() {
  local flowcell
  declare -A seen=()
  declare -a normalized=()

  if [[ ${#FLOWCELLS[@]} -eq 0 ]]; then
    FLOWCELLS=("${DEFAULT_FLOWCELLS[@]}")
  fi

  for flowcell in "${FLOWCELLS[@]}"; do
    [[ "$flowcell" =~ ^[A-Za-z0-9._-]+$ ]] || die "Invalid flowcell name: $flowcell"
    if [[ -z "${seen[$flowcell]:-}" ]]; then
      normalized+=("$flowcell")
      seen["$flowcell"]=1
    fi
  done

  FLOWCELLS=("${normalized[@]}")
}

track_dir() {
  local dir_path="$1"

  if [[ -d "$dir_path" ]]; then
    return
  fi

  mkdir -p "$dir_path"
}

normalize_root_dir() {
  mkdir -p "$ROOT_DIR"
  ROOT_DIR="$(cd "$ROOT_DIR" && pwd -P)"
}

link_path() {
  local source_path="$1"
  local target_path="$2"

  mkdir -p "$(dirname "$source_path")"

  if [[ -e "$source_path" && ! -L "$source_path" ]]; then
    return
  fi

  ln -sfn "$target_path" "$source_path"
}

create_workspace() {
  local suffix
  local flowcell
  local stage_dir

  track_dir "$ROOT_DIR"

  for suffix in "${GLOBAL_DIR_SUFFIXES[@]}"; do
    track_dir "$ROOT_DIR/$suffix"
  done

  track_dir "$ROOT_DIR/$RUNS_ROOT"

  for flowcell in "${FLOWCELLS[@]}"; do
    track_dir "$ROOT_DIR/$RUNS_ROOT/$flowcell"

    for stage_dir in "${FLOWCELL_STAGE_DIRS[@]}"; do
      track_dir "$ROOT_DIR/$RUNS_ROOT/$flowcell/$stage_dir"
    done
  done
}

create_legacy_aliases() {
  local flowcell
  local alias_pair
  local legacy_name
  local standard_name
  local source_path
  local target_path

  for alias_pair in "${LEGACY_ROOT_ALIASES[@]}"; do
    legacy_name="${alias_pair%%:*}"
    standard_name="${alias_pair#*:}"
    source_path="$ROOT_DIR/$legacy_name"
    target_path="$ROOT_DIR/$standard_name"
    link_path "$source_path" "$target_path"
  done

  for flowcell in "${FLOWCELLS[@]}"; do
    source_path="$ROOT_DIR/$flowcell"
    target_path="$ROOT_DIR/$RUNS_ROOT/$flowcell"
    link_path "$source_path" "$target_path"

    for alias_pair in "${LEGACY_FLOWCELL_ALIASES[@]}"; do
      legacy_name="${alias_pair%%:*}"
      standard_name="${alias_pair#*:}"
      source_path="$ROOT_DIR/$RUNS_ROOT/$flowcell/$legacy_name"
      target_path="$ROOT_DIR/$RUNS_ROOT/$flowcell/$standard_name"
      link_path "$source_path" "$target_path"
    done

    source_path="$ROOT_DIR/${flowcell}_p1_filtered_fastq"
    target_path="$ROOT_DIR/$RUNS_ROOT/$flowcell/filtered_reads"
    link_path "$source_path" "$target_path"
  done
}

sync_module4_reference_links() {
  local source_root
  local source_roots=(
    "$ROOT_DIR/results/annotation/corrections/mitogenomes_medaka_rrnS"
    "$ROOT_DIR/results/annotation/mitogenomes_medaka_rrnS"
  )
  local fasta_path
  local fasta_name
  local fasta_stem
  local flowcell
  local barcode
  local target_dir
  local target_path
  declare -A linked=()

  for source_root in "${source_roots[@]}"; do
    [[ -d "$source_root" ]] || continue

    while IFS= read -r -d '' fasta_path; do
      fasta_name="$(basename "$fasta_path")"
      fasta_stem="${fasta_name%.*}"

      if [[ ! "$fasta_stem" =~ ^(F[0-9]+)_barcode([0-9]+)(_.+)?$ ]]; then
        continue
      fi

      flowcell="${BASH_REMATCH[1]}"
      barcode="barcode${BASH_REMATCH[2]}"
      target_dir="$ROOT_DIR/$RUNS_ROOT/$flowcell/references/reoriented_mitogenomes"
      target_path="$target_dir/${barcode}_reoriented.fasta"

      track_dir "$target_dir"

      if [[ -n "${linked[$target_path]:-}" ]]; then
        continue
      fi

      if [[ ! -e "$target_path" || -L "$target_path" ]]; then
        ln -sfn "$fasta_path" "$target_path"
      fi

      linked["$target_path"]=1
    done < <(
      find "$source_root" -type f \
        \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) \
        -print0
    )
  done
}

download_zenodo_data() {
  local zenodo_dir="$ROOT_DIR/downloads/zenodo_fastq"

  track_dir "$zenodo_dir"

  # Add the Zenodo download here when available.
  :
}

main() {
  parse_args "$@"
  [[ -n "$ROOT_DIR" ]] || die "A workspace root is required. Use --root or pass it as the first argument."

  normalize_root_dir
  normalize_flowcells
  create_workspace
  create_legacy_aliases
  sync_module4_reference_links
  download_zenodo_data
}

main "$@"
