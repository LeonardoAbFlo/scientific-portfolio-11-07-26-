source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

set -euo pipefail

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32

for SAMPLE in "$SAMPLE1" "$SAMPLE2"; do
  bcftools stats ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_filtered.vcf.gz \
    > ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_variant_stats.txt
  
  bcftools view -v snps ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_filtered.vcf.gz \
    | bcftools query -f '.\n' \
    | wc -l \
    > ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_snp_count.txt
  
  bcftools view -v indels ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_filtered.vcf.gz \
    | bcftools query -f '.\n' \
    | wc -l \
    > ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_indel_count.txt
done
