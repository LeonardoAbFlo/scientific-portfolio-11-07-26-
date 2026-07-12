source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

set -euo pipefail

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32

for SAMPLE in "$SAMPLE1" "$SAMPLE2"; do
  gatk VariantRecalibrator \
    -R ${REFERENCE} \
    -V ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_snps_recalibrated.vcf.gz \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 ${ROOT_W}/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${ROOT_W}/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -an QD -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL \
    --max-gaussians 4 \
    -O ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_indels.recal \
    --tranches-file ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_indels.tranches

  gatk ApplyVQSR \
    -R ${REFERENCE} \
    -V ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_snps_recalibrated.vcf.gz \
    --recal-file ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_indels.recal \
    --tranches-file ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_indels.tranches \
    -mode INDEL \
    --truth-sensitivity-filter-level 95.0 \
    -O ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_filtered.vcf.gz
done
