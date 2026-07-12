source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32

for SAMPLE in "$SAMPLE1" "$SAMPLE2"; do
  IN_GVCF="${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}.g.vcf.gz"
  OUT_VCF="${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_raw_variants.vcf.gz"
  OUT_DIR="${ROOT_W}/results/var/${SAMPLE}"
  mkdir -p "$OUT_DIR"

  if [[ -s "$OUT_VCF" && ( -s "${OUT_VCF}.tbi" || -s "${OUT_VCF}.csi" ) ]]; then
    echo "[skip] ${SAMPLE}: GenotypeGVCFs already exists"
    continue
  fi

  gatk --java-options "-Xmx16g" GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V "$IN_GVCF" \
    -O "$OUT_VCF"
    
  if [[ ! -s "${OUT_VCF}.tbi" && ! -s "${OUT_VCF}.csi" ]]; then
    tabix -p vcf "$OUT_VCF"
  fi
done
