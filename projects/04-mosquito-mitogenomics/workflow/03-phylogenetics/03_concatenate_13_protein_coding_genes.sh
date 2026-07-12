#!/usr/bin/env bash
set -euo pipefail

# Concatenate the 13 mitochondrial PCGs from the prepared phylogeny dataset.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd -P)}"
PHYLO_ROOT="${PHYLO_ROOT:-$WORKSPACE_ROOT/results/phylogeny}"
REPORT_ROOT="${REPORT_ROOT:-$WORKSPACE_ROOT/results/reports/module3}"
LOG_DIR="${LOG_DIR:-$WORKSPACE_ROOT/logs}"
LOG_PATH="${LOG_DIR}/03_concatenate_13_protein_coding_genes.log"

ENV_NAME="${ENV_NAME:-mito_phylo310}"
RESULTS_DIR="${RESULTS_DIR:-$PHYLO_ROOT/mitofinder_results}"
CONCAT_DIR="${CONCAT_DIR:-$PHYLO_ROOT/concat}"
REPORT="${REPORT:-$REPORT_ROOT/03_concatenate_13_protein_coding_genes.report.tsv}"

GENE_ORDER=(cox1 cox2 cox3 cytb atp6 atp8 nad1 nad2 nad3 nad4 nad4l nad5 nad6)

ALIAS_JSON='{
  "cob":"cytb",
  "nd1":"nad1","nd2":"nad2","nd3":"nad3",
  "nd4":"nad4","nd4l":"nad4l",
  "nd5":"nad5","nd6":"nad6"
}'

mkdir -p "$CONCAT_DIR" "$REPORT_ROOT" "$LOG_DIR"
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

[[ -d "$RESULTS_DIR" ]] || die "Phylogeny dataset directory does not exist: $RESULTS_DIR"

echo "========================================"
echo "[START] 13-PCG concatenation"
echo "Input:"
echo "  $RESULTS_DIR"
echo
echo "Output:"
echo "  $CONCAT_DIR"
echo
echo "Environment:"
echo "  $ENV_NAME"
echo
echo "Report:"
echo "  $REPORT"
echo "========================================"

activate_env "$ENV_NAME"
command -v python3 >/dev/null 2>&1 || die "python3 is not available in the active environment."

AA_OUT="$CONCAT_DIR/13pcg_aa.fasta"
NT_OUT="$CONCAT_DIR/13pcg_nt.fasta"

printf 'sample\taa_source\tnt_source\tstatus\tnote\n' > "$REPORT"

python3 - "$RESULTS_DIR" "$AA_OUT" "$NT_OUT" "$REPORT" "${GENE_ORDER[*]}" "$ALIAS_JSON" <<'PY'
import collections
import glob
import json
import os
import re
import sys
import textwrap
from pathlib import Path

from Bio.Seq import Seq

root, aa_out, nt_out, report_path, order, alias_json = sys.argv[1:7]
order = order.split()
aliases = json.loads(alias_json)


def canon(name: str) -> str:
    name = name.lower()
    name = re.sub(r"\(.*?\)$", "", name)
    name = re.sub(r"_[0-9]+$", "", name)
    return aliases.get(name, name)


def read_fasta(path: str) -> dict[str, str]:
    records: collections.OrderedDict[str, list[str]] = collections.OrderedDict()
    header = None
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip()
            if line.startswith(">"):
                header = canon(line.split("@")[-1])
                records[header] = []
            elif header is not None:
                records[header].append(line)
    return {key: "".join(value) for key, value in records.items()}


def wrap(sequence: str) -> str:
    return "\n".join(textwrap.wrap(sequence, 80))


aa_fastas = sorted(glob.glob(os.path.join(root, "**", "*_final_genes_AA.fasta"), recursive=True))
if not aa_fastas:
    raise SystemExit(f"No *_final_genes_AA.fasta files were found in {root}")

report_rows = []
accepted = 0

with open(aa_out, "w") as aa_handle, open(nt_out, "w") as nt_handle:
    for aa_fasta in aa_fastas:
        sample = os.path.basename(aa_fasta).split("_final")[0]
        nt_fasta = aa_fasta.replace("_AA.fasta", "_NT.fasta")

        if not os.path.isfile(nt_fasta):
            print(f"[WARN] Missing NT FASTA for {sample}")
            report_rows.append((sample, aa_fasta, nt_fasta, "MISSING_NT", "NT FASTA not found"))
            continue

        aa_records = read_fasta(aa_fasta)
        nt_records = read_fasta(nt_fasta)

        aa_concat = []
        nt_concat = []
        status = "OK"
        note = ""

        for gene in order:
            if gene not in aa_records or gene not in nt_records:
                status = "MISSING_GENE"
                note = f"Missing {gene}"
                print(f"[WARN] Missing gene {gene} in {sample}")
                break

            aa_seq = aa_records[gene]
            nt_seq = nt_records[gene]

            if aa_seq.endswith("*"):
                aa_seq = aa_seq[:-1]
                nt_seq = nt_seq[:-3]

            if "*" in aa_seq:
                status = "INTERNAL_STOP"
                note = f"Internal stop in {gene}"
                print(f"[WARN] Internal stop in {sample}:{gene}; sample excluded")
                break

            aa_concat.append(aa_seq)
            nt_concat.append(nt_seq)

        if status == "OK":
            aa_joined = "".join(aa_concat)
            nt_joined = "".join(nt_concat)

            if len(nt_joined) % 3 != 0 or str(Seq(nt_joined).translate(table=5)) != aa_joined:
                status = "AA_NT_MISMATCH"
                note = "Translated NT sequence does not match the AA sequence"
                print(f"[WARN] AA/NT mismatch in {sample}; sample excluded")
            else:
                aa_handle.write(f">{sample}\n{wrap(aa_joined)}\n")
                nt_handle.write(f">{sample}\n{wrap(nt_joined)}\n")
                accepted += 1
                if accepted % 5 == 0 or accepted == len(aa_fastas):
                    print(f"[INFO] Accepted {accepted}/{len(aa_fastas)} samples")

        report_rows.append((sample, aa_fasta, nt_fasta, status, note))

with open(report_path, "a") as report_handle:
    for row in report_rows:
        report_handle.write("\t".join(row) + "\n")

print(f"[INFO] Final matrices written for {accepted}/{len(aa_fastas)} samples")
PY

echo
echo "========================================"
echo "[DONE] 13-PCG concatenation complete"
echo "AA matrix: $AA_OUT"
echo "NT matrix: $NT_OUT"
echo "Report: $REPORT"
echo "Log: $LOG_PATH"
echo "========================================"
