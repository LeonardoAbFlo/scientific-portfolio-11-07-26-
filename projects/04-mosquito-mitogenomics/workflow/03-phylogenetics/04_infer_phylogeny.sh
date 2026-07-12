#!/usr/bin/env bash
set -euo pipefail

# Build the 13-PCG phylogeny from the prepared concatenated matrices.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd -P)}"
PHYLO_ROOT="${PHYLO_ROOT:-$WORKSPACE_ROOT/results/phylogeny}"
REPORT_ROOT="${REPORT_ROOT:-$WORKSPACE_ROOT/results/reports/module3}"
LOG_DIR="${LOG_DIR:-$WORKSPACE_ROOT/logs}"
LOG_PATH="${LOG_DIR}/04_infer_phylogeny.log"

ENV_NAME="${ENV_NAME:-mito_phylo310}"
THREADS="${THREADS:-20}"
MODEL="${MODEL:-GTR+F+I+G4}"
BOOTSTRAP="${BOOTSTRAP:-1000}"
ALRT="${ALRT:-1000}"

CONCAT_DIR="${CONCAT_DIR:-$PHYLO_ROOT/concat}"
TMP_DIR="${TMP_DIR:-$PHYLO_ROOT/tmp}"
TREE_DIR="${TREE_DIR:-$PHYLO_ROOT/tree}"
PREFIX="${PREFIX:-mito_phylo}"
OUT_PREFIX="${OUT_PREFIX:-$TREE_DIR/$PREFIX}"
REPORT="${REPORT:-$REPORT_ROOT/04_infer_phylogeny.report.tsv}"

mkdir -p "$TMP_DIR" "$TREE_DIR" "$REPORT_ROOT" "$LOG_DIR"
exec > >(tee -a "$LOG_PATH") 2>&1

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

activate_env() {
  local env_name="$1"

  if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook -s bash)"
    micromamba activate "$env_name"
  elif command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$env_name"
  else
    die "Could not find micromamba or conda."
  fi
}

resolve_cmd() {
  local cmd_name="$1"
  shift
  local candidate

  if command -v "$cmd_name" >/dev/null 2>&1; then
    command -v "$cmd_name"
    return 0
  fi

  for candidate in "$@"; do
    if [[ -x "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done

  die "Could not find $cmd_name."
}

AA_INPUT="$CONCAT_DIR/13pcg_aa.fasta"
NT_INPUT="$CONCAT_DIR/13pcg_nt.fasta"
AA_ALN="$TMP_DIR/aa_aln.fasta"
NT_ALN="$TMP_DIR/nt_aln.fasta"

[[ -f "$AA_INPUT" ]] || die "AA matrix does not exist: $AA_INPUT"
[[ -f "$NT_INPUT" ]] || die "NT matrix does not exist: $NT_INPUT"

echo "========================================"
echo "[START] 13-PCG phylogeny"
echo "AA matrix:"
echo "  $AA_INPUT"
echo
echo "NT matrix:"
echo "  $NT_INPUT"
echo
echo "Tree prefix:"
echo "  $OUT_PREFIX"
echo
echo "Environment:"
echo "  $ENV_NAME"
echo "========================================"

activate_env "$ENV_NAME"

MAFFT_BIN="${MAFFT_BIN:-$(resolve_cmd mafft "$HOME/anaconda3/envs/mafft/bin/mafft" "$HOME/anaconda3/envs/mitophy/bin/mafft")}"
PAL2NAL_BIN="${PAL2NAL_BIN:-$(resolve_cmd pal2nal.pl "$HOME/anaconda3/envs/pal2nal/bin/pal2nal.pl")}"
IQTREE_BIN="${IQTREE_BIN:-$(resolve_cmd iqtree "$HOME/anaconda3/envs/iqtree/bin/iqtree" "$HOME/anaconda3/envs/mitophy/bin/iqtree")}"

echo "$(timestamp) [INFO] Running MAFFT"
"$MAFFT_BIN" --thread "$THREADS" --auto --anysymbol "$AA_INPUT" > "$AA_ALN"

echo "$(timestamp) [INFO] Running pal2nal"
"$PAL2NAL_BIN" "$AA_ALN" "$NT_INPUT" -codontable 5 -nomismatch -output fasta > "$NT_ALN"

echo "$(timestamp) [INFO] Running IQ-TREE"
"$IQTREE_BIN" \
  -s "$NT_ALN" \
  -m "$MODEL" \
  -B "$BOOTSTRAP" \
  -alrt "$ALRT" \
  -T "$THREADS" \
  -redo \
  --prefix "$OUT_PREFIX"

printf 'aa_alignment\tnt_alignment\ttree_prefix\tstatus\n%s\t%s\t%s\tOK\n' "$AA_ALN" "$NT_ALN" "$OUT_PREFIX" > "$REPORT"

echo
echo "========================================"
echo "[DONE] 13-PCG phylogeny complete"
echo "AA alignment: $AA_ALN"
echo "NT alignment: $NT_ALN"
echo "Tree prefix: $OUT_PREFIX"
echo "Report: $REPORT"
echo "Log: $LOG_PATH"
echo "========================================"
