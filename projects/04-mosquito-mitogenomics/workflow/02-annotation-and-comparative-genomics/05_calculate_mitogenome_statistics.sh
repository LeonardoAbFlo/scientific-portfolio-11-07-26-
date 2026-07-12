#!/usr/bin/env bash
set -euo pipefail

# Summarize mapping and coverage statistics for curated mitogenomes.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd -P)}"
RUNS_ROOT="${RUNS_ROOT:-$WORKSPACE_ROOT/runs}"
ANNOTATION_ROOT="${ANNOTATION_ROOT:-$WORKSPACE_ROOT/results/annotation}"
REPORT_ROOT="${REPORT_ROOT:-$WORKSPACE_ROOT/results/reports/module2}"

THREADS="${THREADS:-30}"
MAP_PRESET="${MAP_PRESET:-map-ont}"
MITO_DIR="${MITO_DIR:-${ANNOTATION_ROOT}/corrections/mitogenomes_medaka_rrnS}"
OUT_TSV="${OUT_TSV:-${REPORT_ROOT}/mitogenome_stats.tsv}"

for cmd in minimap2 samtools seqkit awk find sort; do
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "ERROR: Could not find '$cmd' in PATH." >&2
    exit 1
  }
done

if [[ ! -d "$MITO_DIR" ]]; then
  echo "ERROR: Mitogenome directory does not exist: $MITO_DIR" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT_TSV")"

FLOWCELLS=("$@")

mapfile -t FASTA_FILES < <(find "$MITO_DIR" -type f -name '*.fasta' | sort -V)

if [[ "${#FASTA_FILES[@]}" -eq 0 ]]; then
  echo "ERROR: No FASTA files were found in $MITO_DIR" >&2
  exit 1
fi

printf "Sample\ttotal_reads\trecruited_reads\trecruited_reads_mapped_to_final_mitogenome\tmitogenome_length_bp\tpct_mapped_of_total_reads\tpct_mapped_of_recruited_reads\tmean_depth\n" > "$OUT_TSV"

for fasta in "${FASTA_FILES[@]}"; do
  fasta_base="$(basename "$fasta" _rrnS.fasta)"

  if [[ "$fasta_base" =~ ^(F[0-9]+)_barcode([0-9]+)(_|$) ]]; then
    flowcell="${BASH_REMATCH[1]}"
    barcode="barcode${BASH_REMATCH[2]}"
  elif [[ "$fasta_base" =~ ^(F[0-9]+)_([0-9]+)(-|_|$) ]]; then
    flowcell="${BASH_REMATCH[1]}"
    barcode="barcode${BASH_REMATCH[2]}"
  else
    echo "WARNING: Unrecognized FASTA name, skipping: $fasta" >&2
    continue
  fi

  sample="${flowcell}_${barcode}"

  if [[ "${#FLOWCELLS[@]}" -gt 0 ]]; then
    keep_sample=false
    for selected_flowcell in "${FLOWCELLS[@]}"; do
      if [[ "$flowcell" == "$selected_flowcell" ]]; then
        keep_sample=true
        break
      fi
    done

    if [[ "$keep_sample" != true ]]; then
      continue
    fi
  fi

  raw_fastq="${RUNS_ROOT}/${flowcell}/reads/${barcode}.fastq.gz"
  recruited_fastq="${RUNS_ROOT}/${flowcell}/filtered_reads/${barcode}.fastq.gz"

  if [[ ! -f "$raw_fastq" ]]; then
    echo "WARNING: Missing raw FASTQ for ${sample}: $raw_fastq" >&2
    continue
  fi

  if [[ ! -f "$recruited_fastq" ]]; then
    echo "WARNING: Missing recruited FASTQ for ${sample}: $recruited_fastq" >&2
    continue
  fi

  total_reads="$(seqkit stats -T "$raw_fastq" | awk '
    NR==1 {
      for(i=1;i<=NF;i++) if($i=="num_seqs") c=i
    }
    NR==2 {
      gsub(/,/, "", $c)
      print $c
    }
  ')"

  recruited_reads="$(seqkit stats -T "$recruited_fastq" | awk '
    NR==1 {
      for(i=1;i<=NF;i++) if($i=="num_seqs") c=i
    }
    NR==2 {
      gsub(/,/, "", $c)
      print $c
    }
  ')"

  recruited_reads_mapped_to_final_mitogenome="$(
    minimap2 -ax "$MAP_PRESET" -t "$THREADS" "$fasta" "$recruited_fastq" 2>/dev/null \
      | samtools view -@ "$THREADS" -c -F 2308 -
  )"

  mito_len="$(awk '
    /^>/ {next}
    {
      gsub(/[[:space:]]/, "", $0)
      L += length($0)
    }
    END {print L+0}
  ' "$fasta")"

  pct_mapped_of_total_reads="$(
    awk -v a="$recruited_reads_mapped_to_final_mitogenome" -v t="$total_reads" \
      'BEGIN{if(t==0) print "0.00"; else printf "%.4f", (a/t)*100}'
  )"

  pct_mapped_of_recruited_reads="$(
    awk -v a="$recruited_reads_mapped_to_final_mitogenome" -v t="$recruited_reads" \
      'BEGIN{if(t==0) print "0.00"; else printf "%.2f", (a/t)*100}'
  )"

  mean_depth="$(
    minimap2 -ax "$MAP_PRESET" -t "$THREADS" "$fasta" "$recruited_fastq" 2>/dev/null \
      | samtools view -@ "$THREADS" -b -F 2308 - \
      | samtools sort -@ "$THREADS" -o - \
      | samtools depth -a - \
      | awk '{sum+=$3; n++} END{if(n==0) print "0.00"; else printf "%.2f", sum/n}'
  )"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sample" \
    "$total_reads" \
    "$recruited_reads" \
    "$recruited_reads_mapped_to_final_mitogenome" \
    "$mito_len" \
    "$pct_mapped_of_total_reads" \
    "$pct_mapped_of_recruited_reads" \
    "$mean_depth" >> "$OUT_TSV"

done

echo "TSV written to: $OUT_TSV"
