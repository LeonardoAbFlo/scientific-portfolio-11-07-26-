#!/usr/bin/env bash
# Run Flye across a coverage grid and keep the best assembly for one barcode.
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <FLOWCELL> <BARCODE_NUMBER>"
  echo "Example: $0 F3 92"
  exit 1
fi

FLOWCELL="$1"
BARCODE_NUM="$2"
BARCODE_NAME="barcode${BARCODE_NUM}"

THREADS=30
CONDA_ENV="flye"
CONDA_SH="/path/to/anaconda3/etc/profile.d/conda.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$SCRIPT_DIR}"
RUNS_ROOT="${RUNS_ROOT:-$WORKSPACE_ROOT/runs}"

FASTQ_FILE="${RUNS_ROOT}/${FLOWCELL}/filtered_reads/${BARCODE_NAME}.fastq.gz"
GRID_OUT_BASE_DIR="${WORKSPACE_ROOT}/work/assemblies/${FLOWCELL}/${BARCODE_NAME}"
FINAL_SELECTED_BASE_DIR="${WORKSPACE_ROOT}/results/assemblies/selected"
FINAL_SELECTED_DIR="${FINAL_SELECTED_BASE_DIR}/${FLOWCELL}/${BARCODE_NAME}"

GENOME_SIZE="16k"
MIN_OVERLAP=1000
COVERAGES=(10 15 20 25 30 35 40 45 50 55 60 65 70 75 80)
MIN_LEN=14500
MAX_LEN=17000
TARGET_LEN=16000

LOG_DIR="${WORKSPACE_ROOT}/logs"
mkdir -p "$LOG_DIR"

SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
LOG_PATH="${LOG_DIR}/${SCRIPT_NAME%.sh}_${FLOWCELL}_${BARCODE_NAME}.log"

exec > >(tee -a "$LOG_PATH") 2>&1

echo "[START] Flye coverage-grid selection"
echo "[INFO] Sample: ${FLOWCELL}/${BARCODE_NAME}"
echo "[INFO] Input FASTQ: $FASTQ_FILE"
echo "[INFO] Grid output: $GRID_OUT_BASE_DIR"
echo "[INFO] Final output: $FINAL_SELECTED_DIR"
echo "[INFO] Threads: $THREADS"
echo "[INFO] Log: $LOG_PATH"

if [[ -f "$CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$CONDA_SH"
else
  echo "ERROR: Missing $CONDA_SH"
  exit 1
fi

conda activate "$CONDA_ENV" || {
  echo "ERROR: Cannot activate conda env $CONDA_ENV"
  exit 1
}

command -v flye >/dev/null 2>&1 || {
  echo "ERROR: flye not found in PATH inside env $CONDA_ENV"
  exit 1
}

command -v awk >/dev/null 2>&1 || {
  echo "ERROR: awk not found"
  exit 1
}

command -v seqkit >/dev/null 2>&1 || {
  echo "WARNING: seqkit not found. FASTA stats will be skipped."
}

if [[ ! -f "$FASTQ_FILE" ]]; then
  echo "ERROR: Input FASTQ not found: $FASTQ_FILE"
  exit 1
fi

mkdir -p "$GRID_OUT_BASE_DIR" "$FINAL_SELECTED_BASE_DIR"

SUMMARY="${GRID_OUT_BASE_DIR}/flye_coverage_grid_selection_summary.tsv"
CANDIDATES="${GRID_OUT_BASE_DIR}/flye_coverage_grid_candidates.tsv"
SELECTION_REPORT="${GRID_OUT_BASE_DIR}/flye_selection_report.txt"

printf 'coverage\tstatus\tcontig\tlength\tcoverage_depth\tcirc\tassembly_fasta\tassembly_info\tflye_out_dir\n' > "$SUMMARY"
printf 'coverage\tcontig\tlength\tdistance_from_target\tcoverage_depth\tcirc\tassembly_fasta\tassembly_info\tflye_out_dir\n' > "$CANDIDATES"

for COV in "${COVERAGES[@]}"; do
  echo
  echo "[INFO] Running Flye with --asm-coverage ${COV}"

  OUT_DIR="${GRID_OUT_BASE_DIR}/cov${COV}"
  ASM="${OUT_DIR}/assembly.fasta"
  INFO="${OUT_DIR}/assembly_info.txt"

  if [[ -d "$OUT_DIR" ]]; then
    echo "[INFO] Removing previous output directory: $OUT_DIR"
    rm -rf "$OUT_DIR"
  fi

  mkdir -p "$OUT_DIR"
  START=$(date +%s)

  if flye \
      --nano-hq "$FASTQ_FILE" \
      --out-dir "$OUT_DIR" \
      --threads "$THREADS" \
      --genome-size "$GENOME_SIZE" \
      --asm-coverage "$COV" \
      --min-overlap "$MIN_OVERLAP"; then
    END=$(date +%s)
    DURATION_MIN=$(awk "BEGIN {printf \"%.2f\", ($END - $START)/60}")
    echo "[OK] Flye completed for --asm-coverage ${COV} in ${DURATION_MIN} minutes."
  else
    echo "ERROR: Flye failed for --asm-coverage ${COV}"
    printf '%s\tfailed\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n' "$COV" "$ASM" "$INFO" "$OUT_DIR" >> "$SUMMARY"
    continue
  fi

  if [[ ! -s "$INFO" ]]; then
    echo "[WARN] assembly_info.txt is missing or empty for coverage ${COV}"
    printf '%s\tmissing_assembly_info\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n' "$COV" "$ASM" "$INFO" "$OUT_DIR" >> "$SUMMARY"
    continue
  fi

  if [[ ! -s "$ASM" ]]; then
    echo "[WARN] assembly.fasta is missing or empty for coverage ${COV}"
    printf '%s\tmissing_assembly_fasta\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n' "$COV" "$ASM" "$INFO" "$OUT_DIR" >> "$SUMMARY"
    continue
  fi

  echo "[INFO] assembly_info.txt"
  cat "$INFO"

  awk -v cov_run="$COV" -v asm="$ASM" -v info="$INFO" -v outdir="$OUT_DIR" '
    BEGIN {
      FS=OFS="\t";
      len_col=0; circ_col=0; cov_col=0; name_col=1;
    }
    NR==1 {
      for (i=1; i<=NF; i++) {
        gsub(/^#/, "", $i);
        if ($i == "seq_name") name_col=i;
        if ($i == "length") len_col=i;
        if ($i == "cov." || $i == "cov") cov_col=i;
        if ($i == "circ." || $i == "circ") circ_col=i;
      }
      next;
    }
    len_col > 0 && circ_col > 0 {
      contig=$name_col;
      len=$len_col;
      depth=(cov_col > 0 ? $cov_col : "NA");
      circ=$circ_col;
      print cov_run, "completed", contig, len, depth, circ, asm, info, outdir;
    }
  ' "$INFO" >> "$SUMMARY"

  awk \
    -v cov_run="$COV" \
    -v asm="$ASM" \
    -v info="$INFO" \
    -v outdir="$OUT_DIR" \
    -v min_len="$MIN_LEN" \
    -v max_len="$MAX_LEN" \
    -v target_len="$TARGET_LEN" '
    BEGIN {
      FS=OFS="\t";
      len_col=0; circ_col=0; cov_col=0; name_col=1;
    }
    NR==1 {
      for (i=1; i<=NF; i++) {
        gsub(/^#/, "", $i);
        if ($i == "seq_name") name_col=i;
        if ($i == "length") len_col=i;
        if ($i == "cov." || $i == "cov") cov_col=i;
        if ($i == "circ." || $i == "circ") circ_col=i;
      }
      next;
    }
    len_col > 0 && circ_col > 0 {
      contig=$name_col;
      len=$len_col + 0;
      depth=(cov_col > 0 ? $cov_col : "NA");
      circ=$circ_col;

      if (len >= min_len && len <= max_len && circ == "Y") {
        dist = len - target_len;
        if (dist < 0) dist = -dist;
        print cov_run, contig, len, dist, depth, circ, asm, info, outdir;
      }
    }
  ' "$INFO" >> "$CANDIDATES"
done

echo
echo "[INFO] Selecting the best candidate"

NUM_CANDIDATES=$(awk 'NR > 1 {n++} END {print n+0}' "$CANDIDATES")

if [[ "$NUM_CANDIDATES" -eq 0 ]]; then
  echo "[WARN] No circular contigs found in the expected length range ${MIN_LEN}-${MAX_LEN} bp."

  {
    echo "No selected assembly."
    echo "Reason: no circular contigs found in expected length range."
    echo "Summary: $SUMMARY"
    echo "Candidates: $CANDIDATES"
  } > "$SELECTION_REPORT"

  echo
  echo "========================================"
  echo "[DONE] Flye selection complete"
  echo "Result: none"
  echo "Summary: $SUMMARY"
  echo "Log: $LOG_PATH"
  echo "========================================"
  exit 0
fi

echo "[INFO] Candidate count: $NUM_CANDIDATES"
echo "[INFO] Candidate table"
column -t -s $'\t' "$CANDIDATES" || cat "$CANDIDATES"

SELECTED_LINE=$(
  awk -v target="$TARGET_LEN" '
    BEGIN { FS=OFS="\t" }
    NR==1 { next }
    {
      cov[NR]=$1
      contig[NR]=$2
      len[NR]=$3 + 0
      dist[NR]=$4 + 0
      depth[NR]=$5 + 0
      circ[NR]=$6
      asm[NR]=$7
      info[NR]=$8
      outdir[NR]=$9

      freq[len[NR]]++
      n=NR
    }
    END {
      if (n < 2) exit 1

      max_freq=0
      for (L in freq) {
        if (freq[L] > max_freq) {
          max_freq=freq[L]
        }
      }

      if (max_freq >= 2) {
        chosen_len=""
        chosen_len_dist=""

        for (L in freq) {
          if (freq[L] == max_freq) {
            d = L - target
            if (d < 0) d = -d

            if (chosen_len == "" || d < chosen_len_dist) {
              chosen_len=L
              chosen_len_dist=d
            }
          }
        }

        best_i=""
        best_depth=-1

        for (i=2; i<=n; i++) {
          if (len[i] == chosen_len) {
            if (depth[i] > best_depth) {
              best_depth=depth[i]
              best_i=i
            }
          }
        }

        print "repeated_length_then_highest_depth", cov[best_i], contig[best_i], len[best_i], dist[best_i], depth[best_i], circ[best_i], asm[best_i], info[best_i], outdir[best_i], max_freq
      } else {
        best_i=""
        best_dist=""
        best_depth=-1

        for (i=2; i<=n; i++) {
          if (best_i == "" || dist[i] < best_dist || (dist[i] == best_dist && depth[i] > best_depth)) {
            best_i=i
            best_dist=dist[i]
            best_depth=depth[i]
          }
        }

        print "closest_to_target_then_highest_depth", cov[best_i], contig[best_i], len[best_i], dist[best_i], depth[best_i], circ[best_i], asm[best_i], info[best_i], outdir[best_i], 1
      }
    }
  ' "$CANDIDATES"
)

SELECTION_RULE="$(echo "$SELECTED_LINE" | awk -F'\t' '{print $1}')"
SELECTED_COV="$(echo "$SELECTED_LINE" | awk -F'\t' '{print $2}')"
SELECTED_CONTIG="$(echo "$SELECTED_LINE" | awk -F'\t' '{print $3}')"
SELECTED_LEN="$(echo "$SELECTED_LINE" | awk -F'\t' '{print $4}')"
SELECTED_DIST="$(echo "$SELECTED_LINE" | awk -F'\t' '{print $5}')"
SELECTED_DEPTH="$(echo "$SELECTED_LINE" | awk -F'\t' '{print $6}')"
SELECTED_CIRC="$(echo "$SELECTED_LINE" | awk -F'\t' '{print $7}')"
SELECTED_ASM="$(echo "$SELECTED_LINE" | awk -F'\t' '{print $8}')"
SELECTED_INFO="$(echo "$SELECTED_LINE" | awk -F'\t' '{print $9}')"
SELECTED_OUT_DIR="$(echo "$SELECTED_LINE" | awk -F'\t' '{print $10}')"
SELECTED_LENGTH_FREQ="$(echo "$SELECTED_LINE" | awk -F'\t' '{print $11}')"

echo "[OK] Selected assembly"
echo "[INFO] Rule: $SELECTION_RULE"
echo "[INFO] Coverage: $SELECTED_COV"
echo "[INFO] Contig: $SELECTED_CONTIG"
echo "[INFO] Length: $SELECTED_LEN bp"
echo "[INFO] Distance from target: $SELECTED_DIST bp"
echo "[INFO] Depth: $SELECTED_DEPTH"
echo "[INFO] Circular: $SELECTED_CIRC"
echo "[INFO] Length frequency: $SELECTED_LENGTH_FREQ"
echo "[INFO] Source: $SELECTED_OUT_DIR"

echo
echo "[INFO] Copying selected Flye output to $FINAL_SELECTED_DIR"

if [[ -d "$FINAL_SELECTED_DIR" ]]; then
  echo "[INFO] Removing previous selected directory: $FINAL_SELECTED_DIR"
  rm -rf "$FINAL_SELECTED_DIR"
fi

mkdir -p "$(dirname "$FINAL_SELECTED_DIR")"
cp -a "$SELECTED_OUT_DIR" "$FINAL_SELECTED_DIR"

{
  echo "sample=${FLOWCELL}/${BARCODE_NAME}"
  echo "selection_rule=${SELECTION_RULE}"
  echo "selected_asm_coverage=${SELECTED_COV}"
  echo "selected_contig=${SELECTED_CONTIG}"
  echo "selected_length_bp=${SELECTED_LEN}"
  echo "distance_from_target_bp=${SELECTED_DIST}"
  echo "target_length_bp=${TARGET_LEN}"
  echo "expected_length_range_bp=${MIN_LEN}-${MAX_LEN}"
  echo "flye_coverage_depth=${SELECTED_DEPTH}"
  echo "circular_flag=${SELECTED_CIRC}"
  echo "length_frequency_among_candidates=${SELECTED_LENGTH_FREQ}"
  echo "source_flye_output_dir=${SELECTED_OUT_DIR}"
  echo "source_assembly_fasta=${SELECTED_ASM}"
  echo "source_assembly_info=${SELECTED_INFO}"
  echo "summary_table=${SUMMARY}"
  echo "candidate_table=${CANDIDATES}"
  echo "log=${LOG_PATH}"
} > "${FINAL_SELECTED_DIR}/selection_metadata.txt"

cp "$SUMMARY" "${FINAL_SELECTED_DIR}/flye_coverage_grid_selection_summary.tsv"
cp "$CANDIDATES" "${FINAL_SELECTED_DIR}/flye_coverage_grid_candidates.tsv"
cp "$SELECTION_REPORT" "${FINAL_SELECTED_DIR}/flye_selection_report.txt" 2>/dev/null || true

{
  echo "Flye mitochondrial assembly selection report"
  echo "==========================================="
  echo
  echo "Sample: ${FLOWCELL}/${BARCODE_NAME}"
  echo "Selection rule: ${SELECTION_RULE}"
  echo "Selected --asm-coverage: ${SELECTED_COV}"
  echo "Selected contig: ${SELECTED_CONTIG}"
  echo "Selected length: ${SELECTED_LEN} bp"
  echo "Target length: ${TARGET_LEN} bp"
  echo "Distance from target: ${SELECTED_DIST} bp"
  echo "Flye coverage depth: ${SELECTED_DEPTH}"
  echo "Circular flag: ${SELECTED_CIRC}"
  echo "Length frequency among candidates: ${SELECTED_LENGTH_FREQ}"
  echo
  echo "Selected Flye output copied to:"
  echo "${FINAL_SELECTED_DIR}"
  echo
  echo "Source Flye output:"
  echo "${SELECTED_OUT_DIR}"
  echo
  echo "Summary table:"
  echo "${SUMMARY}"
  echo
  echo "Candidate table:"
  echo "${CANDIDATES}"
} > "${FINAL_SELECTED_DIR}/flye_selection_report.txt"

if command -v seqkit >/dev/null 2>&1; then
  echo
  echo "[INFO] Selected assembly stats"
  seqkit stats "${FINAL_SELECTED_DIR}/assembly.fasta" || true
fi

echo
echo "========================================"
echo "[DONE] Flye selection complete"
echo "Selected assembly: $FINAL_SELECTED_DIR"
echo "Summary: $SUMMARY"
echo "Candidates: $CANDIDATES"
echo "Log: $LOG_PATH"
echo "========================================"
