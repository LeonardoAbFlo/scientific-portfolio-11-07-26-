#!/usr/bin/env bash
set -euo pipefail
GMX=${GMX:-gmx_mpi}
DISC_PS=${DISC_PS:-10000}

for d in rep{1..5}; do
  echo "=== $d ==="
  tpr="$d/md.tpr"
  xtc="$d/md_fit.xtc"
  ndx="$d/index.ndx"
  [[ -f "$tpr" && -f "$xtc" && -f "$ndx" ]] || { echo "[skip] missing files"; continue; }

  printf "RNA\nRNA\n" | $GMX rms -s "$tpr" -f "$xtc" -n "$ndx" -o "$d/rmsd_rna.xvg" -b "$DISC_PS" >/dev/null

  printf "RNA\n" | $GMX gyrate -s "$tpr" -f "$xtc" -n "$ndx" -o "$d/rg_rna.xvg" -b "$DISC_PS" >/dev/null

  printf "RNA\n" | $GMX rmsf -s "$tpr" -f "$xtc" -n "$ndx" -o "$d/rmsf_rna.xvg" -res -b "$DISC_PS" >/dev/null

  echo "[ok] rmsd/rg/rmsf in $d"
done

# Use with: DISC_PS=10000 bash 08_metrics.sh
