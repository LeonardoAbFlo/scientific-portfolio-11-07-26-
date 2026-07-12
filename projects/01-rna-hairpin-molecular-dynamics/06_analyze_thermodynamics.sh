#!/usr/bin/env bash
set -euo pipefail
GMX=${GMX:-gmx_mpi}

usage() {
  echo "Usage:"
  echo "  $0 --list <in.edr>"
  echo "  $0 <in.edr> <out.xvg> <term1> [term2 ...]"
  echo
  echo "Example:"
  echo "  $0 --list em.edr | head"
  echo "  $0 em.edr em_thermo.xvg Potential \"Total\""
}

list_terms() {
  local edr="$1"
  [[ -f "$edr" ]] || { echo "Missing: $edr" >&2; exit 1; }

  printf "0\n" | "$GMX" energy -f "$edr" -o /dev/null 2>&1 \
  | awk '
    function flush() {
      if (idx != "" && name != "") {
        gsub(/^[ \t]+|[ \t]+$/, "", name)
        print idx "\t" name
      }
      idx=""; name=""
    }
    /^[ \t]*[0-9]+[ \t]/ {
      for (i=1; i<=NF; i++) {
        if ($i ~ /^[0-9]+$/) {
          flush()
          idx=$i
          name=""
        } else if (idx != "") {
          name = (name=="" ? $i : name " " $i)
        }
      }
      flush()
    }'
}

find_index() {
  local edr="$1"
  local pat="$2"

  local hits
  hits="$(list_terms "$edr" | awk -v p="$pat" 'BEGIN{IGNORECASE=1} $2 ~ p {print $1 "\t" $2}')"

  local n
  n="$(printf "%s\n" "$hits" | sed '/^$/d' | wc -l | tr -d ' ')"

  if [[ "$n" -eq 0 ]]; then
    echo ""  # caller will warn
    return 0
  fi

  if [[ "$n" -gt 1 ]]; then
    echo "[error] Pattern '$pat' matches multiple terms in $edr:" >&2
    printf "%s\n" "$hits" >&2
    echo "[hint] Use a more specific pattern (e.g., 'Potential Energy' vs 'Potential')." >&2
    exit 2
  fi

  printf "%s\n" "$hits" | awk '{print $1}'
}

if [[ $# -lt 1 ]]; then usage; exit 1; fi

if [[ "${1:-}" == "--list" ]]; then
  [[ $# -eq 2 ]] || { usage; exit 1; }
  list_terms "$2"
  exit 0
fi

if [[ $# -lt 3 ]]; then usage; exit 1; fi

edr="$1"; out="$2"; shift 2
terms=("$@")

idxs=()
for t in "${terms[@]}"; do
  idx="$(find_index "$edr" "$t" || true)"
  if [[ -n "${idx:-}" ]]; then
    idxs+=("$idx")
  else
    echo "[warn] Term not found in $edr (pattern): '$t'" >&2
  fi
done

if [[ "${#idxs[@]}" -eq 0 ]]; then
  echo "[skip] No terms matched. Run: $0 --list $edr" >&2
  exit 3
fi

{ for i in "${idxs[@]}"; do echo "$i"; done; echo 0; } \
  | "$GMX" energy -f "$edr" -o "$out" >/dev/null

echo "[ok] wrote $out"

# Examples:
# ./06_analyze_thermodynamics.sh em.edr em_potential.xvg "Potential"
# ./06_analyze_thermodynamics.sh nvt.edr nvt_temperature.xvg "Temperature"
# ./06_analyze_thermodynamics.sh npt.edr npt_pressure.xvg "Pressure"
# ./06_analyze_thermodynamics.sh npt.edr npt_density.xvg "Density"
# ./06_analyze_thermodynamics.sh npt.edr npt_volume.xvg "Volume"
