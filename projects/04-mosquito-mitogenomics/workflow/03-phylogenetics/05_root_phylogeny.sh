#!/usr/bin/env bash
set -euo pipefail

# Re-root the phylogeny with the configured outgroup label.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd -P)}"
PHYLO_ROOT="${PHYLO_ROOT:-$WORKSPACE_ROOT/results/phylogeny}"
REPORT_ROOT="${REPORT_ROOT:-$WORKSPACE_ROOT/results/reports/module3}"
LOG_DIR="${LOG_DIR:-$WORKSPACE_ROOT/logs}"
LOG_PATH="${LOG_DIR}/05_root_phylogeny.log"

ENV_NAME="${ENV_NAME:-mito_phylo310}"
TREE_DIR="${TREE_DIR:-$PHYLO_ROOT/tree}"
TREE_PREFIX="${TREE_PREFIX:-$TREE_DIR/mito_phylo}"
IN_TREE="${IN_TREE:-${TREE_PREFIX}.treefile}"
OUT_TREE="${OUT_TREE:-${TREE_PREFIX}_rooted.tree}"
OUT_ID="${OUT_ID:-outgroup}"
REPORT="${REPORT:-$REPORT_ROOT/05_root_phylogeny.report.tsv}"

mkdir -p "$TREE_DIR" "$REPORT_ROOT" "$LOG_DIR"
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

[[ -f "$IN_TREE" ]] || die "Input tree does not exist: $IN_TREE"

echo "========================================"
echo "[START] Tree re-rooting"
echo "Input tree:"
echo "  $IN_TREE"
echo
echo "Outgroup:"
echo "  $OUT_ID"
echo
echo "Output tree:"
echo "  $OUT_TREE"
echo "========================================"

activate_env "$ENV_NAME"

if command -v nw_reroot >/dev/null 2>&1; then
  nw_reroot "$IN_TREE" "$OUT_ID" > "$OUT_TREE"
else
  python3 - "$OUT_ID" "$IN_TREE" "$OUT_TREE" <<'PY'
import sys
from Bio import Phylo

out_id, in_tree, out_tree = sys.argv[1:4]
tree = Phylo.read(in_tree, "newick")
target = next((clade for clade in tree.find_clades() if clade.name == out_id), None)
if target is None:
    raise SystemExit(f"Outgroup was not found in the tree: {out_id}")
tree.root_with_outgroup(target)
Phylo.write(tree, out_tree, "newick")
PY
fi

printf 'input_tree\toutgroup\toutput_tree\tstatus\n%s\t%s\t%s\tOK\n' "$IN_TREE" "$OUT_ID" "$OUT_TREE" > "$REPORT"

echo
echo "========================================"
echo "[DONE] Tree re-rooting complete"
echo "Output tree: $OUT_TREE"
echo "Report: $REPORT"
echo "Log: $LOG_PATH"
echo "========================================"
