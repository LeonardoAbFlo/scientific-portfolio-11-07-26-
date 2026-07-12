#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=scripts/lib/common.sh
source "$SCRIPT_DIR/scripts/lib/common.sh"

usage() { printf 'Usage: run_methylome_pipeline.sh --config FILE --samples TSV [--force]\n'; }
config=""; samples=""; force=0
while (($#)); do
  case "$1" in
    --config) config="$2"; shift 2 ;;
    --samples) samples="$2"; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done
[[ -n "$config" && -f "$config" && -n "$samples" && -f "$samples" ]] || { usage >&2; exit 2; }
# shellcheck disable=SC1090
source "$config"
: "${POD5_DIR:?Set POD5_DIR in the config}"
: "${WORK_DIR:?Set WORK_DIR in the config}"
: "${DORADO_MODEL:?Set DORADO_MODEL in the config}"
: "${DORADO_KIT:?Set DORADO_KIT in the config}"
THREADS="${THREADS:-16}"
DORADO_DEVICE="${DORADO_DEVICE:-auto}"
MIN_QSCORE="${MIN_QSCORE:-10}"
MIN_COVERAGE="${MIN_COVERAGE:-20}"
CANDIDATE_FRACTION="${CANDIDATE_FRACTION:-0.20}"
WEAK_FRACTION="${WEAK_FRACTION:-0.05}"

force_arg=()
((force)) && force_arg+=(--force)
"$SCRIPT_DIR/scripts/01_basecall_modified_bases.sh" --pod5 "$POD5_DIR" --output "$WORK_DIR/basecalls" --model "$DORADO_MODEL" --kit "$DORADO_KIT" --device "$DORADO_DEVICE" --min-qscore "$MIN_QSCORE" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/02_demultiplex_reads.sh" --basecalls "$WORK_DIR/basecalls" --output "$WORK_DIR/demultiplexed" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/03_convert_bam_to_fastq.sh" --bams "$WORK_DIR/demultiplexed" --output "$WORK_DIR/reads" --threads "$THREADS" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/04_align_reads_to_mitogenomes.sh" --reads "$WORK_DIR/reads" --samples "$samples" --output "$WORK_DIR/alignments" --threads "$THREADS" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/05_call_methylated_sites.sh" --alignments "$WORK_DIR/alignments" --output "$WORK_DIR/pileups" --threads "$THREADS" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/06_discover_methylation_motifs.sh" --pileups "$WORK_DIR/pileups" --samples "$samples" --output "$WORK_DIR/motifs" --threads "$THREADS" "${force_arg[@]}"
"$SCRIPT_DIR/scripts/07_summarize_methylation_results.py" --samples "$samples" --pileups "$WORK_DIR/pileups" --motifs "$WORK_DIR/motifs" --output "$WORK_DIR/reports" --min-coverage "$MIN_COVERAGE" --candidate-fraction "$CANDIDATE_FRACTION" --weak-fraction "$WEAK_FRACTION"
log "Methylome workflow complete: $WORK_DIR/reports"
