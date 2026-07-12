source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32

for SAMPLE in "$SAMPLE1" "$SAMPLE2"; do
  IN_BAM="${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}_recalibrated.bam"
  QC_DIR="${ROOT_W}/results/qc/${SAMPLE}"
  mkdir -p "$QC_DIR"

  gatk CollectAlignmentSummaryMetrics \
    -R "${REFERENCE}" \
    -I "$IN_BAM" \
    -O "${QC_DIR}/${SAMPLE}_alignment_summary.txt"

  gatk CollectInsertSizeMetrics \
    -I "$IN_BAM" \
    -O "${QC_DIR}/${SAMPLE}_insert_size_metrics.txt" \
    -H "${QC_DIR}/${SAMPLE}_insert_size_histogram.pdf"
done
