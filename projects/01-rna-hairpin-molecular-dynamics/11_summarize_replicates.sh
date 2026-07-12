#!/usr/bin/env bash
set -euo pipefail

avg2() { awk '($1!~/#/ && $1!~/@/){s+=$2;n++} END{if(n>0) printf "%.4f", s/n; else printf "NA"}' "$1"; }
frac_gt0() { awk '($1!~/#/ && $1!~/@/){n++; if($2>0) k++} END{if(n>0) printf "%.3f", k/n; else printf "NA"}' "$1"; }

echo -e "rep\tavg_HB\tfrac(HB>0)\tavg_contacts\tfrac(contacts>0)\tavg_end2end"
for d in rep{1..5}; do
  hb="$d/hb_stem.xvg"
  ct="$d/contacts_stem.xvg"
  ee="$d/end2end_com.xvg"
  echo -e "$d\t$(avg2 "$hb")\t$(frac_gt0 "$hb")\t$(avg2 "$ct")\t$(frac_gt0 "$ct")\t$(avg2 "$ee")"
done
