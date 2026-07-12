source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

set -euo pipefail

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32

for SAMPLE in "$SAMPLE1" "$SAMPLE2"; do
  snpEff ann \
  -Xmx32g \
  -stats ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_annotation_stats.html \
  GRCh38.105 \
  ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_filtered.vcf.gz \
  > ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_annotated.vcf

  gatk VariantsToTable \
    -V ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_annotated.vcf \
    -F CHROM -F POS -F REF -F ALT -F QUAL -F ANN \
    -O ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_variants_table.tsv
done
