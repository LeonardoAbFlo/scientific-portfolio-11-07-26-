#!/usr/bin/env bash
set -euo pipefail
GMX=${GMX:-gmx_mpi}
DISC_PS=${DISC_PS:-10000}
pairs=("1 17" "2 16" "3 15" "4 14" "5 13")

for d in rep{1..5}; do
  tpr="$d/md.tpr"; xtc="$d/md_fit.xtc"
  [[ -f "$tpr" && -f "$xtc" ]] || continue
  for p in "${pairs[@]}"; do
    a=${p% *}; b=${p#* }
    $GMX distance -s "$tpr" -f "$xtc" -select "com of resid $a plus com of resid $b" \
      -oall "$d/bp_${a}_${b}.xvg" -b "$DISC_PS" >/dev/null
  done
done

# Use with: DISC_PS=10000 bash 13-hairpin-distance.sh
