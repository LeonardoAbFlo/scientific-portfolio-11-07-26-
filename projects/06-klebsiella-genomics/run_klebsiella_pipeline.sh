#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=scripts/lib/common.sh
source "$SCRIPT_DIR/scripts/lib/common.sh"

usage() {
  printf 'Usage: run_klebsiella_pipeline.sh --config FILE [--force]\n'
}

config=""; force=0
while (($#)); do
  case "$1" in
    --config) config="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done
[[ -n "$config" && -f "$config" ]] || { usage >&2; exit 2; }

# shellcheck disable=SC1090
source "$config"
: "${READS_DIR:?Set READS_DIR in the config}"
: "${WORK_DIR:?Set WORK_DIR in the config}"
: "${METADATA_TSV:?Set METADATA_TSV in the config}"
THREADS="${THREADS:-16}"
FLYE_MODE="${FLYE_MODE:-nano-hq}"

force_arg=()
((force)) && force_arg+=(--force)

"$SCRIPT_DIR/scripts/01_assemble_reads.sh" --reads "$READS_DIR" --output "$WORK_DIR/assemblies" --threads "$THREADS" --mode "$FLYE_MODE" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/02_split_contigs.py" --assemblies "$WORK_DIR/assemblies" --output "$WORK_DIR/contigs" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/03_profile_assemblies.sh" --assemblies "$WORK_DIR/assemblies" --output "$WORK_DIR/initial_profiles" --threads "$THREADS" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/04_classify_contigs.sh" --contigs "$WORK_DIR/contigs" --output "$WORK_DIR/contig_qc" --threads "$THREADS" --gtdb-data "${GTDBTK_DATA_PATH:?Set GTDBTK_DATA_PATH}" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/05_build_taxonomy_table.py" --manifest "$WORK_DIR/contigs/contig_manifest.tsv" --mobsuite "$WORK_DIR/initial_profiles/mobsuite" --gtdbtk "$WORK_DIR/contig_qc/gtdbtk" --output "$WORK_DIR/tables/taxonomy.tsv"

filter_args=()
[[ "${EXCLUDE_CROSS_GENUS_PLASMIDS:-0}" == 1 ]] && filter_args+=(--exclude-cross-genus-plasmids)
"$SCRIPT_DIR/scripts/06_filter_contigs.py" --metadata "$METADATA_TSV" --taxonomy "$WORK_DIR/tables/taxonomy.tsv" --contigs "$WORK_DIR/contigs" --output "$WORK_DIR/filtered_contigs" "${filter_args[@]}" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/07_merge_contigs.py" --contigs "$WORK_DIR/filtered_contigs" --output "$WORK_DIR/filtered_assemblies" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/08_profile_filtered_assemblies.sh" --assemblies "$WORK_DIR/filtered_assemblies" --output "$WORK_DIR/filtered_profiles" --threads "$THREADS" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/09_assess_read_quality.sh" --reads "$READS_DIR" --output "$WORK_DIR/read_qc" --threads "$THREADS"

log "Computational workflow complete. Review filter_decisions.tsv before generating reports."
