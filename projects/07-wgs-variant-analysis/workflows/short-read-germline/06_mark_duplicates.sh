source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32

for SAMPLE in "$SAMPLE1" "$SAMPLE2"; do
  IN_BAM="${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}.bam"
  OUT_BAM="${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}_marked_duplicates.bam"
  METRICS="${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}_duplicate_metrics.txt"

  if [[ -s "$OUT_BAM" && -s "${OUT_BAM}.bai" ]]; then
    echo "[skip] ${SAMPLE}: MarkDuplicates already done"
    continue
  fi

  gatk MarkDuplicates \
    -I "$IN_BAM" \
    -O "$OUT_BAM" \
    -M "$METRICS" \
    --CREATE_INDEX true
done
