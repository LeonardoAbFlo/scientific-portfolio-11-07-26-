source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32
DBSNP="${ROOT_W}/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz"

for SAMPLE in "$SAMPLE1" "$SAMPLE2"; do
  IN_BAM="${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}_recalibrated.bam"
  OUTDIR="${ROOT_W}/results/var/${SAMPLE}"
  OUTGVCF="${OUTDIR}/${SAMPLE}.g.vcf.gz"
  mkdir -p "$OUTDIR"
  
  if [[ -s "$OUTGVCF" && -s "${OUTGVCF}.tbi" ]]; then
    echo "[skip] ${SAMPLE}: gVCF already exists"
    continue
  fi

  gatk --java-options "-Xmx16g" HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "$IN_BAM" \
    -O "$OUTGVCF" \
    -ERC GVCF \
    --dbsnp "$DBSNP" \
    --native-pair-hmm-threads "${THREADS}"
done
