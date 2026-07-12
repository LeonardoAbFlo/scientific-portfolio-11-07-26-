#!/usr/bin/env bash
set -euo pipefail
GMX=${GMX:-gmx_mpi}

for d in rep{1..5}; do
  echo "=== $d ==="
  gro="$d/npt.gro"
  ndx="$d/index.ndx"
  out="$d/index_hairpin.ndx"
  [[ -f "$gro" && -f "$ndx" ]] || { echo "[skip] missing $gro or $ndx"; continue; }

  ngrp=$(grep -c '^\[' "$ndx")
  STEM5_ID=$ngrp
  STEM3_ID=$((ngrp+1))
  LOOP_ID=$((ngrp+2))
  RES1_ID=$((ngrp+3))
  RES17_ID=$((ngrp+4))

  {
    echo "r 1-5"
    echo "r 13-17"
    echo "r 6-12"
    echo "r 1"
    echo "r 17"
    echo "name $STEM5_ID STEM5"
    echo "name $STEM3_ID STEM3"
    echo "name $LOOP_ID LOOP"
    echo "name $RES1_ID RES1"
    echo "name $RES17_ID RES17"
    echo "q"
  } | $GMX make_ndx -f "$gro" -n "$ndx" -o "$out" >/dev/null

  echo "[ok] wrote $out"
done
