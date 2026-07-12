#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob dotglob

# Assemble the combined phylogeny dataset from module 2 results and phylogeny references.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd -P)}"
ANNOTATION_ROOT="${ANNOTATION_ROOT:-$WORKSPACE_ROOT/results/annotation}"
PHYLO_ROOT="${PHYLO_ROOT:-$WORKSPACE_ROOT/results/phylogeny}"
REPORT_ROOT="${REPORT_ROOT:-$WORKSPACE_ROOT/results/reports/module3}"
LOG_DIR="${LOG_DIR:-$WORKSPACE_ROOT/logs}"
LOG_PATH="${LOG_DIR}/02_prepare_phylogeny_dataset.log"

SAMPLE_RESULTS_ROOT="${SAMPLE_RESULTS_ROOT:-$ANNOTATION_ROOT/corrections/medaka_mitofinder_rrnS}"
REFERENCE_RESULTS_ROOT="${REFERENCE_RESULTS_ROOT:-$PHYLO_ROOT/reference_mitofinder}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$PHYLO_ROOT/mitofinder_results}"
MODE="${MODE:-link}"
FORCE="${FORCE:-0}"

REPORT="${REPORT:-$REPORT_ROOT/02_prepare_phylogeny_dataset.report.tsv}"

usage() {
  cat <<'EOF'
Usage:
  02_prepare_phylogeny_dataset.sh

Key variables:
  SAMPLE_RESULTS_ROOT=/path/to/results/annotation/corrections/medaka_mitofinder_rrnS
  REFERENCE_RESULTS_ROOT=/path/to/results/phylogeny/reference_mitofinder
  OUTPUT_ROOT=/path/to/results/phylogeny/mitofinder_results
  MODE=link
  FORCE=1
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

mkdir -p "$OUTPUT_ROOT" "$REPORT_ROOT" "$LOG_DIR"
exec > >(tee -a "$LOG_PATH") 2>&1

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

record_report() {
  printf "%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$5" >> "$REPORT"
}

stage_item() {
  local source_group="$1"
  local source_path="$2"
  local target_path="$3"

  if [[ ! -d "$source_path" ]]; then
    record_report "$source_group" "$source_path" "$target_path" "missing" "Source directory was not found."
    return 1
  fi

  if [[ -e "$target_path" ]]; then
    if [[ "$FORCE" == "1" ]]; then
      rm -rf -- "$target_path"
    else
      record_report "$source_group" "$source_path" "$target_path" "skipped" "Target already exists."
      return 0
    fi
  fi

  mkdir -p "$(dirname "$target_path")"

  case "$MODE" in
    link)
      ln -s "$source_path" "$target_path"
      record_report "$source_group" "$source_path" "$target_path" "linked" "Directory linked into the phylogeny dataset."
      ;;
    copy)
      cp -a "$source_path" "$target_path"
      record_report "$source_group" "$source_path" "$target_path" "copied" "Directory copied into the phylogeny dataset."
      ;;
    *)
      die "Unsupported MODE: $MODE"
      ;;
  esac
}

[[ -d "$SAMPLE_RESULTS_ROOT" ]] || die "Sample result directory not found: $SAMPLE_RESULTS_ROOT"
[[ -d "$REFERENCE_RESULTS_ROOT" ]] || die "Reference result directory not found: $REFERENCE_RESULTS_ROOT"

printf "source_group\tsource_path\ttarget_path\tstatus\tnote\n" > "$REPORT"

echo "[$(timestamp)] Preparing combined phylogeny dataset"
echo "[$(timestamp)] Sample results:    $SAMPLE_RESULTS_ROOT"
echo "[$(timestamp)] Reference results: $REFERENCE_RESULTS_ROOT"
echo "[$(timestamp)] Output root:       $OUTPUT_ROOT"
echo "[$(timestamp)] Mode:              $MODE"

sample_count=0
reference_count=0

for genus_dir in "$SAMPLE_RESULTS_ROOT"/*; do
  [[ -d "$genus_dir" ]] || continue
  genus_name="$(basename "$genus_dir")"
  mkdir -p "$OUTPUT_ROOT/$genus_name"

  for sample_dir in "$genus_dir"/*; do
    [[ -d "$sample_dir" ]] || continue
    sample_name="$(basename "$sample_dir")"
    stage_item "sample" "$sample_dir" "$OUTPUT_ROOT/$genus_name/$sample_name"
    sample_count=$((sample_count + 1))
  done
done

for reference_dir in "$REFERENCE_RESULTS_ROOT"/*; do
  [[ -d "$reference_dir" ]] || continue
  reference_name="$(basename "$reference_dir")"

  case "$reference_name" in
    anopheles|culex|psorophora)
      record_report "reference" "$reference_dir" "$OUTPUT_ROOT/$reference_name" "skipped" "Reserved genus name; sample groups are staged from module 2."
      continue
      ;;
  esac

  stage_item "reference" "$reference_dir" "$OUTPUT_ROOT/$reference_name"
  reference_count=$((reference_count + 1))
done

echo "[$(timestamp)] Samples processed:    $sample_count"
echo "[$(timestamp)] References processed: $reference_count"
echo "[$(timestamp)] Report: $REPORT"
echo "[DONE] Phylogeny dataset prepared"
