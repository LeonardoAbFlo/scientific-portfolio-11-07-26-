#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Run MitoFinder on manually corrected rrnS FASTA files.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd -P)}"
ANNOTATION_ROOT="${ANNOTATION_ROOT:-$WORKSPACE_ROOT/results/annotation}"
REPORT_ROOT="${REPORT_ROOT:-$WORKSPACE_ROOT/results/reports/module2}"

IN_ROOT="${IN_ROOT:-$ANNOTATION_ROOT/corrections/mitogenomes_medaka_rrnS}"
OUT_ROOT="${OUT_ROOT:-$ANNOTATION_ROOT/corrections/medaka_mitofinder_rrnS}"

REF_DIR="${REF_DIR:-$ANNOTATION_ROOT/ref_gb}"
REF_GB="${REF_GB:-$REF_DIR/culicidae_mt_refseq.gb}"
REF_GB_ORIGINAL="${REF_GB_ORIGINAL:-}"

GENCODE="${GENCODE:-5}"
THREADS="${THREADS:-32}"
MAXMEM="${MAXMEM:-16}"
CONDA_ENV="${CONDA_ENV:-sif_env}"

REPORT="${REPORT:-$REPORT_ROOT/04_apply_manual_annotation_corrections.report.tsv}"

mkdir -p "$REF_DIR" "$OUT_ROOT" "$(dirname "$REPORT")"

prepare_reference() {
  if [[ -f "$REF_GB" ]]; then
    return 0
  fi

  if [[ -n "$REF_GB_ORIGINAL" && -f "$REF_GB_ORIGINAL" ]]; then
    echo "[INFO] Copying GenBank reference into the workspace"
    echo "  from: $REF_GB_ORIGINAL"
    echo "  to:   $REF_GB"
    cp "$REF_GB_ORIGINAL" "$REF_GB"
    return 0
  fi

  echo "ERROR: GenBank reference not found."
  echo "Checked:"
  echo "  $REF_GB"
  if [[ -n "$REF_GB_ORIGINAL" ]]; then
    echo "  $REF_GB_ORIGINAL"
  fi
  exit 1
}

activate_mitofinder_env() {
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV"
}

prepare_reference

if [[ ! -d "$IN_ROOT" ]]; then
  echo "ERROR: Input directory does not exist:"
  echo "  $IN_ROOT"
  exit 1
fi

activate_mitofinder_env

if ! command -v mitofinder >/dev/null 2>&1; then
  echo "ERROR: mitofinder is not available in the current environment."
  echo "Expected environment: $CONDA_ENV"
  exit 1
fi

cat > "$REPORT" <<'EOF'
genus	sample	input_fasta	output_dir	seqid	status
EOF

echo "========================================"
echo "[START] MitoFinder after rrnS corrections"
echo "Input:"
echo "  $IN_ROOT"
echo
echo "Output:"
echo "  $OUT_ROOT"
echo
echo "Reference:"
echo "  $REF_GB"
echo
echo "Report:"
echo "  $REPORT"
echo "========================================"
echo

processed=0
ok=0
failed=0
missing_fasta=0

for GENUS_DIR in "$IN_ROOT"/*/; do
  GENUS="$(basename "$GENUS_DIR")"

  for SAMPLE_DIR in "$GENUS_DIR"*/; do
    SAMPLE_NAME="$(basename "$SAMPLE_DIR")"

    mapfile -t FASTAS < <(
      find "$SAMPLE_DIR" -maxdepth 1 -type f \
        \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) \
        | sort
    )

    if ((${#FASTAS[@]} == 0)); then
      echo "----------------------------------------"
      echo "Sample: $GENUS/$SAMPLE_NAME"
      echo "  ! No FASTA files were found in:"
      echo "    $SAMPLE_DIR"
      echo -e "${GENUS}\t${SAMPLE_NAME}\tNA\tNA\tNA\tMISSING_FASTA" >> "$REPORT"
      ((missing_fasta+=1))
      echo
      continue
    fi

    if ((${#FASTAS[@]} > 1)); then
      echo "----------------------------------------"
      echo "Sample: $GENUS/$SAMPLE_NAME"
      echo "  ! More than one FASTA was found. Using the first one:"
      printf '    %s\n' "${FASTAS[@]}"
      echo
    fi

    INPUT_FASTA="${FASTAS[0]}"

    SAMPLE_SAFE="$(echo "$SAMPLE_NAME" | sed -E 's/[^A-Za-z0-9_]+/_/g; s/_+$//')"
    FASTA_BASE="$(basename "$INPUT_FASTA")"
    FASTA_STEM="${FASTA_BASE%.*}"

    SEQID="$(echo "$FASTA_STEM" | sed -E 's/[^A-Za-z0-9_]+/_/g; s/_+$//')"

    if [[ -z "$SEQID" ]]; then
      SEQID="${SAMPLE_SAFE}_rrnS"
    fi

    OUT_SAMPLE_DIR="$OUT_ROOT/$GENUS/$SAMPLE_NAME"
    mkdir -p "$OUT_SAMPLE_DIR"

    LOG_FILE="$OUT_SAMPLE_DIR/${SEQID}.mitofinder.log"

    echo "----------------------------------------"
    echo "Genus:    $GENUS"
    echo "Sample:   $SAMPLE_NAME"
    echo "FASTA:    $INPUT_FASTA"
    echo "SEQID:    $SEQID"
    echo "Output:   $OUT_SAMPLE_DIR"
    echo "Log:      $LOG_FILE"
    echo

    if (
      cd "$OUT_SAMPLE_DIR"

      mitofinder \
        --seqid "$SEQID" \
        --assembly "$INPUT_FASTA" \
        --refseq "$REF_GB" \
        --organism "$GENCODE" \
        --processors "$THREADS" \
        --max-memory "$MAXMEM" \
        --override \
        --blast-size 30 \
        --min-contig-size 1000
    ) > "$LOG_FILE" 2>&1; then

      echo "  [OK] MitoFinder completed: $GENUS/$SAMPLE_NAME"
      echo -e "${GENUS}\t${SAMPLE_NAME}\t${INPUT_FASTA}\t${OUT_SAMPLE_DIR}\t${SEQID}\tOK" >> "$REPORT"
      ((ok+=1))

      COI_FASTA="$OUT_SAMPLE_DIR/COI_${SEQID}.fna"

      mapfile -t CDS_FILES < <(
        find "$OUT_SAMPLE_DIR" -type f -path "*/CDS/*.fna" 2>/dev/null | sort
      )

      if ((${#CDS_FILES[@]} > 0)); then
        awk '
          BEGIN { IGNORECASE=1 }
          /^>/ {
            keep = ($0 ~ /(COI|COX1|CO1|cytochrome c oxidase subunit I|cytochrome oxidase subunit I)/)
          }
          keep { print }
        ' "${CDS_FILES[@]}" > "${COI_FASTA}.tmp"

        if [[ -s "${COI_FASTA}.tmp" ]]; then
          mv "${COI_FASTA}.tmp" "$COI_FASTA"
          echo "  [OK] COI/COX1 written to:"
          echo "    $COI_FASTA"
        else
          rm -f "${COI_FASTA}.tmp"
          echo "  ! COI/COX1 was not found in this run"
        fi
      else
        echo "  ! No CDS/*.fna files are available yet"
      fi

    else
      echo "  ! ERROR: MitoFinder failed for $GENUS/$SAMPLE_NAME"
      echo "  Review the log:"
      echo "    $LOG_FILE"
      echo -e "${GENUS}\t${SAMPLE_NAME}\t${INPUT_FASTA}\t${OUT_SAMPLE_DIR}\t${SEQID}\tFAILED" >> "$REPORT"
      ((failed+=1))
    fi

    ((processed+=1))
    echo
  done
done

echo "========================================"
echo "[DONE] MitoFinder after rrnS corrections complete"
echo "Samples processed:    $processed"
echo "Completed:            $ok"
echo "Failed:               $failed"
echo "Missing FASTA:        $missing_fasta"
echo
echo "Results:"
echo "  $OUT_ROOT"
echo
echo "Report:"
echo "  $REPORT"
echo "========================================"
