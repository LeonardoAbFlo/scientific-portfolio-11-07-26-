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
    -V ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_raw_variants.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${ROOT_W}/reference/hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 ${ROOT_W}/reference/1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${ROOT_W}/reference/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${ROOT_W}/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_snps.recal \
    --tranches-file ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_snps.tranches

  gatk ApplyVQSR \
    -R ${REFERENCE} \
    -V ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_raw_variants.vcf.gz \
    --recal-file ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_snps.recal \
    --tranches-file ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_snps.tranches \
    -mode SNP \
    --truth-sensitivity-filter-level 99.0 \
    -O ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_snps_recalibrated.vcf.gz
done
