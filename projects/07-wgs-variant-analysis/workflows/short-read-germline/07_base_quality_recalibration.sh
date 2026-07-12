source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32

DBSNP="${ROOT_W}/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
MILLS="${ROOT_W}/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

for SAMPLE in "$SAMPLE1" "$SAMPLE2"; do
  IN_BAM="${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}_marked_duplicates.bam"
  RECAL_DIR="${ROOT_W}/results/recal/${SAMPLE}"
  TABLE="${RECAL_DIR}/${SAMPLE}_recal_data.table"
  OUT_BAM="${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}_recalibrated.bam"

  mkdir -p "$RECAL_DIR"

  if [[ ! -s "$TABLE" ]]; then
    gatk BaseRecalibrator \
      -I "$IN_BAM" \
      -R "${REFERENCE}" \
      --known-sites "$DBSNP" \
      --known-sites "$MILLS" \
      -O "$TABLE"
  else
    echo "[skip] ${SAMPLE}: recal table already exists"
  fi

  if [[ -s "$OUT_BAM" && -s "${OUT_BAM}.bai" ]]; then
    echo "[skip] ${SAMPLE}: recalibrated BAM already exists"
  else
    gatk ApplyBQSR \
      -I "$IN_BAM" \
      -R "${REFERENCE}" \
      --bqsr-recal-file "$TABLE" \
      -O "$OUT_BAM"
    samtools index -@ "${THREADS}" "$OUT_BAM"
  fi
done
