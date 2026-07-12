#!/usr/bin/env bash
set -Eeuo pipefail

ROOT_W="${ROOT_W:?Set ROOT_W to the workflow root}"
REFERENCE="${REFERENCE:-$ROOT_W/reference/Homo_sapiens_assembly38.fasta}"
SAMPLES="${SAMPLES:-SRR7890824 SRR7890827}"
THREADS="${THREADS:-32}"

for SAMPLE in $SAMPLES; do
  mkdir -p "${ROOT_W}/results/aligned/${SAMPLE}"

  bwa mem \
    -t "$THREADS" \
    -M \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}_lib" \
    "$REFERENCE" \
    "$ROOT_W/results/trimmed/${SAMPLE}/${SAMPLE}_1_val_1.fq.gz" \
    "$ROOT_W/results/trimmed/${SAMPLE}/${SAMPLE}_2_val_2.fq.gz" \
  | samtools sort -@ "$THREADS" -o "$ROOT_W/results/aligned/${SAMPLE}/${SAMPLE}.bam" -

  samtools index -@ "$THREADS" "$ROOT_W/results/aligned/${SAMPLE}/${SAMPLE}.bam"
done
