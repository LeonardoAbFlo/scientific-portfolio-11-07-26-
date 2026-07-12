#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Annotate phylogeny reference genomes and outgroup FASTAs with MitoFinder.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd -P)}"
ANNOTATION_ROOT="${ANNOTATION_ROOT:-$WORKSPACE_ROOT/results/annotation}"
PHYLO_ROOT="${PHYLO_ROOT:-$WORKSPACE_ROOT/results/phylogeny}"
RESOURCE_ROOT="${RESOURCE_ROOT:-$WORKSPACE_ROOT/resources/phylogeny}"
REPORT_ROOT="${REPORT_ROOT:-$WORKSPACE_ROOT/results/reports/module3}"
LOG_DIR="${LOG_DIR:-$WORKSPACE_ROOT/logs}"
LOG_PATH="${LOG_DIR}/01_annotate_phylogeny_references.log"

REFERENCE_FASTA_ROOT="${REFERENCE_FASTA_ROOT:-$RESOURCE_ROOT/reference_fastas}"
OUTGROUP_FASTA_ROOT="${OUTGROUP_FASTA_ROOT:-$RESOURCE_ROOT/outgroup_fastas}"
WORK_ROOT="${WORK_ROOT:-$WORKSPACE_ROOT/tmp/module3_reference_fastas}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$PHYLO_ROOT/reference_mitofinder}"

REF_DIR="${REF_DIR:-$ANNOTATION_ROOT/ref_gb}"
REF_GB="${REF_GB:-$REF_DIR/culicidae_mt_refseq.gb}"
REF_GB_ORIGINAL="${REF_GB_ORIGINAL:-}"

CONDA_ENV="${CONDA_ENV:-sif_env}"
GENCODE="${GENCODE:-5}"
THREADS="${THREADS:-32}"
MAXMEM="${MAXMEM:-16}"
BLAST_SIZE="${BLAST_SIZE:-30}"
MIN_CONTIG_SIZE="${MIN_CONTIG_SIZE:-1000}"
FORCE="${FORCE:-0}"
ADJUST_DIRECTION="${ADJUST_DIRECTION:-1}"

REPORT="${REPORT:-$REPORT_ROOT/01_annotate_phylogeny_references.report.tsv}"

usage() {
  cat <<'EOF'
Usage:
  01_annotate_phylogeny_references.sh

Key variables:
  REFERENCE_FASTA_ROOT=/path/to/reference_fastas
  OUTGROUP_FASTA_ROOT=/path/to/outgroup_fastas
  OUTPUT_ROOT=/path/to/results/phylogeny/reference_mitofinder
  REF_GB=/path/to/culicidae_mt_refseq.gb
  CONDA_ENV=sif_env
  FORCE=1
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

mkdir -p "$OUTPUT_ROOT" "$WORK_ROOT" "$REPORT_ROOT" "$LOG_DIR" "$REF_DIR"
exec > >(tee -a "$LOG_PATH") 2>&1

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
  log "ERROR: $*"
  exit 1
}

prepare_reference() {
  if [[ -f "$REF_GB" ]]; then
    return 0
  fi

  if [[ -n "$REF_GB_ORIGINAL" && -f "$REF_GB_ORIGINAL" ]]; then
    log "Copying GenBank reference into the workspace"
    log "  from: $REF_GB_ORIGINAL"
    log "  to:   $REF_GB"
    cp "$REF_GB_ORIGINAL" "$REF_GB"
    return 0
  fi

  die "GenBank reference not found: $REF_GB"
}

activate_mitofinder_env() {
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV"
}

split_fasta_file() {
  local input_path="$1"
  local target_root="$2"

  mkdir -p "$target_root"

  python3 - "$input_path" "$target_root" <<'PY'
import re
import sys
from pathlib import Path

input_path = Path(sys.argv[1])
target_root = Path(sys.argv[2])

records = []
header = None
seq_parts = []

with input_path.open() as handle:
    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(re.sub(r"\s+", "", line).upper())

if header is not None:
    records.append((header, "".join(seq_parts)))

if not records:
    raise SystemExit(f"No FASTA records were found in {input_path}")

for header, sequence in records:
    seqid = header.split()[0]
    seqid = re.sub(r"[^A-Za-z0-9._-]+", "_", seqid).strip("._-")
    if not seqid:
        raise SystemExit(f"Could not derive a sequence ID from header '{header}' in {input_path}")
    output_path = target_root / f"{seqid}.fasta"
    if output_path.exists():
        raise SystemExit(f"Duplicate sequence ID detected while splitting FASTA inputs: {seqid}")
    with output_path.open("w") as handle:
        handle.write(f">{seqid}\n")
        for index in range(0, len(sequence), 80):
            handle.write(sequence[index:index + 80] + "\n")
PY
}

prepare_input_fastas() {
  local input_root="$1"
  local category="$2"
  local split_root="$WORK_ROOT/$category"
  local found_any=0
  local input_path

  mkdir -p "$split_root"

  if [[ ! -d "$input_root" ]]; then
    log "INFO: Input FASTA directory not found, skipping: $input_root"
    return 0
  fi

  while IFS= read -r -d '' input_path; do
    found_any=1
    split_fasta_file "$input_path" "$split_root"
  done < <(find "$input_root" -maxdepth 1 -type f \( -iname '*.fa' -o -iname '*.fasta' -o -iname '*.fna' \) -print0 | sort -z)

  if [[ "$found_any" -eq 0 ]]; then
    log "INFO: No FASTA files were found in $input_root"
  fi
}

is_current() {
  local output="$1"
  shift
  local input

  [[ -s "$output" ]] || return 1
  for input in "$@"; do
    [[ "$output" -nt "$input" ]] || return 1
  done
}

command -v python3 >/dev/null 2>&1 || die "python3 is required to split input FASTA files."

prepare_reference
prepare_input_fastas "$REFERENCE_FASTA_ROOT" reference
prepare_input_fastas "$OUTGROUP_FASTA_ROOT" outgroup

mapfile -t INPUT_FASTAS < <(find "$WORK_ROOT" -type f -name '*.fasta' | sort)

if [[ "${#INPUT_FASTAS[@]}" -eq 0 ]]; then
  die "No FASTA files are available for phylogeny references or outgroup."
fi

activate_mitofinder_env

if ! command -v mitofinder >/dev/null 2>&1; then
  die "mitofinder is not available in the current environment: $CONDA_ENV"
fi

printf 'category\tseqid\tinput_fasta\toutput_dir\tstatus\n' > "$REPORT"

log "Workspace root: $WORKSPACE_ROOT"
log "Reference FASTA root: $REFERENCE_FASTA_ROOT"
log "Outgroup FASTA root: $OUTGROUP_FASTA_ROOT"
log "Output root: $OUTPUT_ROOT"
log "Reference GB: $REF_GB"
log "Report: $REPORT"

processed=0
skipped_current=0
failed=0

for INPUT_FASTA in "${INPUT_FASTAS[@]}"; do
  CATEGORY="$(basename "$(dirname "$INPUT_FASTA")")"
  SEQID="$(basename "$INPUT_FASTA" .fasta)"
  OUT_DIR="$OUTPUT_ROOT/$SEQID"
  FINAL_AA="$OUT_DIR/$SEQID/${SEQID}_MitoFinder_mitfi_Final_Results/${SEQID}_final_genes_AA.fasta"

  mkdir -p "$OUT_DIR"

  if [[ "$FORCE" != "1" ]] && is_current "$FINAL_AA" "$INPUT_FASTA" "$REF_GB"; then
    log "Skipping current output for $SEQID"
    printf '%s\t%s\t%s\t%s\tskipped_current\n' "$CATEGORY" "$SEQID" "$INPUT_FASTA" "$OUT_DIR" >> "$REPORT"
    ((skipped_current += 1))
    continue
  fi

  log "Running MitoFinder for $SEQID"

  MITOFINDER_CMD=(
    mitofinder
    --seqid "$SEQID"
    --assembly "$INPUT_FASTA"
    -o "$OUT_DIR"
    --refseq "$REF_GB"
    --organism "$GENCODE"
    --processors "$THREADS"
    --max-memory "$MAXMEM"
    --override
    --blast-size "$BLAST_SIZE"
    --min-contig-size "$MIN_CONTIG_SIZE"
  )

  if [[ "$ADJUST_DIRECTION" == "1" ]]; then
    MITOFINDER_CMD+=(--adjust-direction)
  fi

  if "${MITOFINDER_CMD[@]}"; then
    printf '%s\t%s\t%s\t%s\tOK\n' "$CATEGORY" "$SEQID" "$INPUT_FASTA" "$OUT_DIR" >> "$REPORT"
    ((processed += 1))
  else
    printf '%s\t%s\t%s\t%s\tFAILED\n' "$CATEGORY" "$SEQID" "$INPUT_FASTA" "$OUT_DIR" >> "$REPORT"
    ((failed += 1))
  fi
done

echo
echo "========================================"
echo "[DONE] Phylogeny reference annotation complete"
echo "Processed: $processed"
echo "Skipped current: $skipped_current"
echo "Failed: $failed"
echo "Output: $OUTPUT_ROOT"
echo "Report: $REPORT"
echo "Log: $LOG_PATH"
echo "========================================"
