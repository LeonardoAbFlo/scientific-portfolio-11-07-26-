#!/usr/bin/env bash
set -euo pipefail
GMX=${GMX:-gmx_mpi}
DISC_PS=${DISC_PS:-10000}

for d in rep{1..5}; do
  echo "=== $d ==="
  tpr="$d/md.tpr"
  xtc="$d/md_fit.xtc"
  ndx="$d/index_hairpin.ndx"
  [[ -f "$tpr" && -f "$xtc" && -f "$ndx" ]] || { echo "[skip] missing files"; continue; }

  printf "STEM5\nSTEM3\n" | $GMX hbond -s "$tpr" -f "$xtc" -n "$ndx" -num "$d/hb_stem.xvg" -b "$DISC_PS" >/dev/null

  printf "STEM5\nSTEM3\n" | $GMX mindist -s "$tpr" -f "$xtc" -n "$ndx" -d 0.35 \
    -od "$d/mindist_stem.xvg" -on "$d/contacts_stem.xvg" -b "$DISC_PS" >/dev/null
    
  $GMX distance -s "$tpr" -f "$xtc" -n "$ndx" \
    -select 'com of group "RES1" plus com of group "RES17"' \
    -oall "$d/end2end_com.xvg" -b "$DISC_PS" >/dev/null

  echo "[ok] hb/contacts/end2end in $d"
done

# Use with: DISC_PS=10000 bash 10-hairpin_metrics.sh
