#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Organize Medaka consensuses by taxonomy and run MitoFinder.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd -P)}"
RUNS_ROOT="${RUNS_ROOT:-$WORKSPACE_ROOT/runs}"
ANNOTATION_ROOT="${ANNOTATION_ROOT:-$WORKSPACE_ROOT/results/annotation}"

DATA_ROOT="${DATA_ROOT:-$ANNOTATION_ROOT/data}"
MEDAKA_ROOT="${MEDAKA_ROOT:-$WORKSPACE_ROOT/results/polishing/dorado_medaka}"
OUT_MEDAKA="${OUT_MEDAKA:-$ANNOTATION_ROOT/mitogenomes_medaka}"
OUT_MITOFINDER="${OUT_MITOFINDER:-$ANNOTATION_ROOT/medaka_mitofinder}"
REF_DIR="${REF_DIR:-$ANNOTATION_ROOT/ref_gb}"
REF_GB="${REF_GB:-$REF_DIR/culicidae_mt_refseq.gb}"
REF_GB_ORIGINAL="${REF_GB_ORIGINAL:-}"

GENCODE="${GENCODE:-5}"
THREADS="${THREADS:-32}"
MAXMEM="${MAXMEM:-16}"
RUN_MITOFINDER="${RUN_MITOFINDER:-1}"
COPY_MODE="${COPY_MODE:-copy}"
CONDA_ENV="${CONDA_ENV:-sif_env}"

prepare_reference() {
  mkdir -p "$REF_DIR"

  if [[ -f "$REF_GB" ]]; then
    return 0
  fi

  if [[ -n "$REF_GB_ORIGINAL" && -f "$REF_GB_ORIGINAL" ]]; then
    echo "[INFO] Copying GenBank reference into the workspace"
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

if [[ ! -d "$DATA_ROOT" ]]; then
  echo "ERROR: DATA_ROOT does not exist: $DATA_ROOT"
  exit 1
fi

if [[ ! -d "$MEDAKA_ROOT" ]]; then
  echo "ERROR: MEDAKA_ROOT does not exist: $MEDAKA_ROOT"
  exit 1
fi

if [[ ! -f "$REF_GB" ]]; then
  echo "ERROR: REF_GB does not exist: $REF_GB"
  exit 1
fi

mkdir -p "$OUT_MEDAKA" "$OUT_MITOFINDER"

if [[ "$RUN_MITOFINDER" == "1" ]]; then
  activate_mitofinder_env
fi

echo "========================================"
echo "[START] Module 2 MitoFinder setup"
echo "Workspace root: $WORKSPACE_ROOT"
echo "Runs root:      $RUNS_ROOT"
echo "Data root:      $DATA_ROOT"
echo "Medaka root:    $MEDAKA_ROOT"
echo "Medaka output:  $OUT_MEDAKA"
echo "MitoFinder out: $OUT_MITOFINDER"
echo "Reference GB:   $REF_GB"
echo "Copy mode:      $COPY_MODE"
echo "Run MitoFinder: $RUN_MITOFINDER"
echo "========================================"
echo

processed=0
missing_medaka=0
missing_consensus=0
failed_mitofinder=0

for GENUS_DIR in "$DATA_ROOT"/*/; do
  GENUS="$(basename "$GENUS_DIR")"

  if [[ "$GENUS" == "ref_gb" ]]; then
    continue
  fi

  for SAMPLE_DIR in "$GENUS_DIR"*/; do
    SAMPLE_NAME="$(basename "$SAMPLE_DIR")"

    if [[ "$SAMPLE_NAME" =~ ^(F[0-9]+)_([0-9]+)-(.+)$ ]]; then
      FLOWCELL="${BASH_REMATCH[1]}"
      BARCODE_NUM="${BASH_REMATCH[2]}"
      BARCODE_DIR="barcode${BARCODE_NUM}"
    else
      echo "[WARN] Could not parse sample directory: $SAMPLE_NAME"
      continue
    fi

    MEDAKA_BARCODE_DIR="$MEDAKA_ROOT/$FLOWCELL/$BARCODE_DIR"
    MEDAKA_SRC="$MEDAKA_BARCODE_DIR/medaka"
    MEDAKA_DEST="$OUT_MEDAKA/$GENUS/$SAMPLE_NAME"

    echo "----------------------------------------"
    echo "Sample:          $SAMPLE_NAME"
    echo "Genus:           $GENUS"
    echo "Flowcell:        $FLOWCELL"
    echo "Barcode:         $BARCODE_DIR"
    echo "Medaka sample:   $MEDAKA_BARCODE_DIR"
    echo "Medaka source:   $MEDAKA_SRC"
    echo "Medaka target:   $MEDAKA_DEST"

    if [[ ! -d "$MEDAKA_BARCODE_DIR" ]]; then
      echo "  ! Missing polished sample directory for $FLOWCELL/$BARCODE_DIR"
      echo "    Expected: $MEDAKA_BARCODE_DIR"
      ((missing_medaka+=1))
      echo
      continue
    fi

    if [[ ! -d "$MEDAKA_SRC" ]]; then
      echo "  ! Missing medaka/ directory for $FLOWCELL/$BARCODE_DIR"
      echo "    Expected: $MEDAKA_SRC"
      ((missing_medaka+=1))
      echo
      continue
    fi

    mkdir -p "$MEDAKA_DEST"

    if [[ "$COPY_MODE" == "copy" ]]; then
      cp -a "$MEDAKA_SRC"/. "$MEDAKA_DEST"/
    elif [[ "$COPY_MODE" == "symlink" ]]; then
      for f in "$MEDAKA_SRC"/*; do
        ln -sf "$f" "$MEDAKA_DEST/$(basename "$f")"
      done
    else
      echo "ERROR: COPY_MODE must be 'copy' or 'symlink'"
      exit 1
    fi

    CONSENSUS="$MEDAKA_DEST/consensus.fasta"

    if [[ ! -f "$CONSENSUS" ]]; then
      POLISHED_FASTA="$MEDAKA_BARCODE_DIR/${FLOWCELL}_${BARCODE_DIR}.dorado.medaka.consensus.fasta"

      if [[ -f "$POLISHED_FASTA" ]]; then
        echo "  ! consensus.fasta was not copied from medaka/"
        echo "  -> Falling back to:"
        echo "     $POLISHED_FASTA"
        cp -f "$POLISHED_FASTA" "$CONSENSUS"
      else
        echo "  ! Missing consensus.fasta in:"
        echo "    $MEDAKA_DEST"
        echo "  ! Missing fallback polished FASTA:"
        echo "    $POLISHED_FASTA"
        ((missing_consensus+=1))
        echo
        continue
      fi
    fi

    echo "  [OK] Medaka files organized"

    if [[ "$RUN_MITOFINDER" == "1" ]]; then
      MF_OUTDIR="$OUT_MITOFINDER/$GENUS/$SAMPLE_NAME"
      mkdir -p "$MF_OUTDIR"

      SAMPLE_SAFE="$(echo "$SAMPLE_NAME" | sed -E 's/[^A-Za-z0-9_]+/_/g; s/_+$//')"

      echo "  -> Running MitoFinder"
      echo "     Output:   $MF_OUTDIR"
      echo "     SeqID:    $SAMPLE_SAFE"
      echo "     Assembly: $CONSENSUS"

      if ! (
        cd "$MF_OUTDIR"

        mitofinder \
          --seqid "$SAMPLE_SAFE" \
          --assembly "$CONSENSUS" \
          --refseq "$REF_GB" \
          --organism "$GENCODE" \
          --processors "$THREADS" \
          --max-memory "$MAXMEM" \
          --override \
          --blast-size 30 \
          --min-contig-size 1000
      ); then
        echo "  ! ERROR: MitoFinder failed for $SAMPLE_NAME"
        ((failed_mitofinder+=1))
        echo
        continue
      fi

      COI_FASTA="$MF_OUTDIR/COI_${SAMPLE_SAFE}.fna"

      mapfile -t CDS_FILES < <(find "$MF_OUTDIR" -type f -path "*/CDS/*.fna" 2>/dev/null | sort)

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
          echo "  [OK] COI/COX1 written to: $COI_FASTA"
        else
          rm -f "${COI_FASTA}.tmp"
          echo "  ! COI/COX1/CO1 was not found in the MitoFinder CDS output"
        fi
      else
        echo "  ! No CDS/*.fna files were found in:"
        echo "    $MF_OUTDIR"
      fi
    fi

    ((processed+=1))
    echo
  done
done

echo "========================================"
echo "[DONE] Module 2 MitoFinder setup complete"
echo "Samples processed: $processed"
echo "Missing Medaka dirs: $missing_medaka"
echo "Missing consensus FASTA: $missing_consensus"
echo "MitoFinder failures: $failed_mitofinder"
echo
echo "Organized Medaka output:"
echo "  $OUT_MEDAKA"
echo
echo "MitoFinder output:"
echo "  $OUT_MITOFINDER"
echo "========================================"
