#!/usr/bin/env bash
set -euo pipefail
GMX=${GMX:-gmx_mpi}
DISC_PS=${DISC_PS:-10000}

for d in rep{1..5}; do
  echo "=== $d ==="
  tpr="$d/md.tpr"
  xtc="$d/md.xtc"
  ndx="$d/index.ndx"
  [[ -f "$tpr" && -f "$xtc" && -f "$ndx" ]] || { echo "[skip] missing files"; continue; }

  printf "RNA\nSystem\n" | $GMX trjconv -s "$tpr" -f "$xtc" -n "$ndx" \
    -o "$d/md_nopbc.xtc" -pbc mol -center -ur compact >/dev/null

  printf "RNA\nRNA\n" | $GMX trjconv -s "$tpr" -f "$d/md_nopbc.xtc" -n "$ndx" \
    -o "$d/md_fit.xtc" -fit rot+trans >/dev/null

  echo "[ok] $d/md_fit.xtc ready"
done

# Run with: DISC_PS=10000 bash 07-prep_trajectories.sh
