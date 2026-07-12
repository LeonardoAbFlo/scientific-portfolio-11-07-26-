#!/usr/bin/env bash
set -euo pipefail
avg2() { awk '($1!~/#/ && $1!~/@/){s+=$2;n++} END{if(n>0) printf "%.4f", s/n; else printf "NA"}' "$1"; }
echo -e "rep\tavg_RMSD\tavg_Rg"
for d in rep{1..5}; do
  echo -e "$d\t$(avg2 "$d/rmsd_rna.xvg")\t$(avg2 "$d/rg_rna.xvg")"
done
