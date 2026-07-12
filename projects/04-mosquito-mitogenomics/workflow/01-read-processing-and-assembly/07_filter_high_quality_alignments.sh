#!/usr/bin/env bash
# Map filtered reads back to selected Flye assemblies and export MAPQ subsets.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$SCRIPT_DIR}"
RUNS_ROOT="${RUNS_ROOT:-$WORKSPACE_ROOT/runs}"

SELECTED_ROOT="${SELECTED_ROOT:-$WORKSPACE_ROOT/results/assemblies/selected}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$WORKSPACE_ROOT/results/mapping}"
LOG_DIR="${LOG_DIR:-$WORKSPACE_ROOT/logs}"
LOG_PATH="${LOG_DIR}/07_filter_high_quality_alignments.log"
THREADS=30
MAP_PRESET="map-ont"
CONDA_ENV="minimap2"
CONDA_SH="/path/to/anaconda3/etc/profile.d/conda.sh"
DRY_RUN=0

declare -a MAPQ_THRESHOLDS=(60)
declare -a FILTERS=()

usage() {
  cat <<'EOF'
Usage:
  07_filter_high_quality_alignments.sh [options] [F1|F3|F4|F3/barcode92|F3 92 ...]

Options:
  --threads N        Thread count. Default: 30
  --out-dir DIR      Output directory. Default: results/mapping
  --mapq LIST        Comma-separated MAPQ thresholds. Default: 60
  --dry-run          Print commands without writing files
  -h, --help         Show this help message

Without filters, the script processes every assembly under results/assemblies/selected.
Passing a flowcell followed by a number, for example "F3 90", is interpreted as "F3/barcode90".
EOF
}

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --threads)
        [[ $# -ge 2 ]] || die "Missing value for --threads"
        THREADS="$2"
        shift 2
        ;;
      --out-dir)
        [[ $# -ge 2 ]] || die "Missing value for --out-dir"
        OUTPUT_ROOT="$2"
        shift 2
        ;;
      --mapq)
        [[ $# -ge 2 ]] || die "Missing value for --mapq"
        IFS=',' read -r -a MAPQ_THRESHOLDS <<< "$2"
        shift 2
        ;;
      --dry-run)
        DRY_RUN=1
        shift
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      --)
        shift
        while [[ $# -gt 0 ]]; do
          FILTERS+=("$1")
          shift
        done
        ;;
      -*)
        die "Unknown option: $1"
        ;;
      *)
        FILTERS+=("$1")
        shift
        ;;
    esac
  done
}

normalize_filters() {
  local idx=0
  local current next
  local -a normalized=()

  while (( idx < ${#FILTERS[@]} )); do
    current="${FILTERS[idx]}"

    if (( idx + 1 < ${#FILTERS[@]} )); then
      next="${FILTERS[idx + 1]}"
      if [[ "$current" =~ ^[Ff][0-9]+$ && "$next" =~ ^[0-9]+$ ]]; then
        normalized+=("${current}/barcode${next}")
        ((idx+=2))
        continue
      fi
    fi

    normalized+=("$current")
    ((idx+=1))
  done

  FILTERS=("${normalized[@]}")
}

normalize_thresholds() {
  local threshold
  local -a cleaned=()

  for threshold in "${MAPQ_THRESHOLDS[@]}"; do
    [[ "$threshold" =~ ^[0-9]+$ ]] || die "Invalid MAPQ threshold: $threshold"
    cleaned+=("$threshold")
  done

  mapfile -t MAPQ_THRESHOLDS < <(printf '%s\n' "${cleaned[@]}" | awk '!seen[$0]++' | sort -n)
  [[ ${#MAPQ_THRESHOLDS[@]} -gt 0 ]] || die "No valid MAPQ thresholds were provided"
}

maybe_activate_conda() {
  if command -v minimap2 >/dev/null 2>&1 && command -v samtools >/dev/null 2>&1; then
    return
  fi

  if [[ -f "$CONDA_SH" ]]; then
    # shellcheck disable=SC1090
    source "$CONDA_SH"
    if command -v conda >/dev/null 2>&1; then
      conda activate "$CONDA_ENV" >/dev/null 2>&1 || true
    fi
  fi
}

require_tools() {
  maybe_activate_conda

  local cmd
  for cmd in minimap2 samtools awk find sort; do
    command -v "$cmd" >/dev/null 2>&1 || die "Could not find '$cmd' in PATH"
  done
}

matches_filters() {
  local sample_id="$1"
  local flowcell="$2"
  local barcode="$3"
  local filter

  if [[ ${#FILTERS[@]} -eq 0 ]]; then
    return 0
  fi

  for filter in "${FILTERS[@]}"; do
    if [[ "$sample_id" == "$filter" || "$flowcell" == "$filter" || "$barcode" == "$filter" ]]; then
      return 0
    fi
  done

  return 1
}

write_mapq_distribution() {
  local bam="$1"
  local min_mapq="$2"
  local out_tsv="$3"

  if (( DRY_RUN )); then
    log "[DRY-RUN] MAPQ distribution >=${min_mapq} -> ${out_tsv}"
    return
  fi

  {
    printf 'mapq\treads\n'
    samtools view "$bam" \
      | awk -v min_mapq="$min_mapq" '$5 >= min_mapq {counts[$5]++} END {for (q in counts) print q "\t" counts[q]}' \
      | sort -k1,1n
  } > "$out_tsv"
}

collect_depth_metrics() {
  local bam="$1"

  samtools depth -a "$bam" | awk '
    BEGIN {sum=0; total=0; covered=0}
    {sum += $3; total++; if ($3 > 0) covered++}
    END {
      if (total == 0) {
        printf "0\t0\t0\t0\n"
      } else {
        printf "%.2f\t%.6f\t%d\t%d\n", sum / total, covered / total, covered, total
      }
    }
  '
}

append_summary_row() {
  local summary_tsv="$1"
  local sample_id="$2"
  local flowcell="$3"
  local barcode="$4"
  local min_mapq="$5"
  local bam="$6"
  local dist_tsv="$7"

  local mapped_reads mean_depth breadth covered_bases reference_bases

  mapped_reads="$(samtools view -@ "$THREADS" -c "$bam")"
  read -r mean_depth breadth covered_bases reference_bases < <(collect_depth_metrics "$bam")

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample_id" \
    "$flowcell" \
    "$barcode" \
    "$min_mapq" \
    "$mapped_reads" \
    "$mean_depth" \
    "$breadth" \
    "$covered_bases" \
    "$reference_bases" \
    "$dist_tsv" >> "$summary_tsv"
}

main() {
  parse_args "$@"
  normalize_filters
  normalize_thresholds
  require_tools
  mkdir -p "$LOG_DIR"
  exec > >(tee -a "$LOG_PATH") 2>&1

  [[ -d "$SELECTED_ROOT" ]] || die "Selected assemblies directory does not exist: $SELECTED_ROOT"

  local manifest_tsv="${OUTPUT_ROOT}/mapping_manifest.tsv"
  local summary_tsv="${OUTPUT_ROOT}/mapq_summary.tsv"
  local assembly sample_dir flowcell barcode sample_id fastq
  local raw_bam raw_dist raw_log filtered_bam filtered_dist threshold
  local processed=0 skipped_missing_fastq=0 skipped_filters=0

  if (( DRY_RUN )); then
    log "[INFO] Dry-run mode enabled"
  else
    mkdir -p "$OUTPUT_ROOT"
    printf 'sample_id\tflowcell\tbarcode\tfastq\tassembly\tout_dir\n' > "$manifest_tsv"
    printf 'sample_id\tflowcell\tbarcode\tmin_mapq\tmapped_reads\tmean_depth\tbreadth\tcovered_bases\treference_bases\tdistribution_tsv\n' > "$summary_tsv"
  fi

  mapfile -t assemblies < <(
    find "$SELECTED_ROOT" -mindepth 3 -maxdepth 3 -type f -path '*/barcode*/assembly.fasta' | sort -V
  )

  [[ ${#assemblies[@]} -gt 0 ]] || die "No assembly.fasta files were found in $SELECTED_ROOT"

  for assembly in "${assemblies[@]}"; do
    barcode="$(basename "$(dirname "$assembly")")"
    flowcell="$(basename "$(dirname "$(dirname "$assembly")")")"
    sample_id="${flowcell}/${barcode}"

    if ! matches_filters "$sample_id" "$flowcell" "$barcode"; then
      ((skipped_filters+=1))
      continue
    fi

    fastq="${RUNS_ROOT}/${flowcell}/filtered_reads/${barcode}.fastq.gz"
    if [[ ! -f "$fastq" ]]; then
      log "[WARN] Missing FASTQ for ${sample_id}: ${fastq}"
      ((skipped_missing_fastq+=1))
      continue
    fi

    sample_dir="${OUTPUT_ROOT}/${flowcell}/${barcode}"
    raw_bam="${sample_dir}/${barcode}.mapped.sorted.bam"
    raw_dist="${sample_dir}/${barcode}.mapq_distribution.tsv"
    raw_log="${sample_dir}/${barcode}.minimap2.stderr.log"

    log "[INFO] Processing ${sample_id}"
    log "  assembly: ${assembly}"
    log "  fastq: ${fastq}"
    log "  output: ${sample_dir}"

    if (( DRY_RUN )); then
      log "[DRY-RUN] mkdir -p ${sample_dir}"
      log "[DRY-RUN] minimap2 -ax ${MAP_PRESET} -t ${THREADS} ${assembly} ${fastq} | samtools view -b -F 4 | samtools sort -o ${raw_bam}"
      log "[DRY-RUN] samtools index ${raw_bam}"
      write_mapq_distribution "$raw_bam" 0 "$raw_dist"
      for threshold in "${MAPQ_THRESHOLDS[@]}"; do
        filtered_bam="${sample_dir}/${barcode}.mapq${threshold}.sorted.bam"
        filtered_dist="${sample_dir}/${barcode}.mapq_distribution_q${threshold}.tsv"
        log "[DRY-RUN] samtools view -h -q ${threshold} ${raw_bam} | samtools sort -o ${filtered_bam}"
        log "[DRY-RUN] samtools index ${filtered_bam}"
        write_mapq_distribution "$filtered_bam" "$threshold" "$filtered_dist"
      done
      ((processed+=1))
      continue
    fi

    mkdir -p "$sample_dir"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$sample_id" "$flowcell" "$barcode" "$fastq" "$assembly" "$sample_dir" >> "$manifest_tsv"

    minimap2 -ax "$MAP_PRESET" -t "$THREADS" "$assembly" "$fastq" 2> "$raw_log" \
      | samtools view -@ "$THREADS" -b -F 4 - \
      | samtools sort -@ "$THREADS" -o "$raw_bam" -

    samtools index -@ "$THREADS" "$raw_bam"
    write_mapq_distribution "$raw_bam" 0 "$raw_dist"
    append_summary_row "$summary_tsv" "$sample_id" "$flowcell" "$barcode" 0 "$raw_bam" "$raw_dist"

    for threshold in "${MAPQ_THRESHOLDS[@]}"; do
      filtered_bam="${sample_dir}/${barcode}.mapq${threshold}.sorted.bam"
      filtered_dist="${sample_dir}/${barcode}.mapq_distribution_q${threshold}.tsv"

      samtools view -@ "$THREADS" -h -q "$threshold" "$raw_bam" \
        | samtools sort -@ "$THREADS" -o "$filtered_bam" -

      if [[ "$(samtools view -@ "$THREADS" -c "$filtered_bam")" -gt 0 ]]; then
        samtools index -@ "$THREADS" "$filtered_bam"
      fi

      write_mapq_distribution "$filtered_bam" "$threshold" "$filtered_dist"
      append_summary_row "$summary_tsv" "$sample_id" "$flowcell" "$barcode" "$threshold" "$filtered_bam" "$filtered_dist"
    done

    ((processed+=1))
  done

  echo
  echo "========================================"
  echo "[DONE] MAPQ mapping complete"
  echo "Processed: ${processed}"
  echo "Skipped by filter: ${skipped_filters}"
  echo "Skipped missing FASTQ: ${skipped_missing_fastq}"
  echo "Output: ${OUTPUT_ROOT}"
  if (( ! DRY_RUN )); then
    echo "Manifest: ${manifest_tsv}"
    echo "Summary: ${summary_tsv}"
  fi
  echo "Log: ${LOG_PATH}"
  echo "========================================"
}

main "$@"
