source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

set -euo pipefail

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32

for SAMPLE in "$SAMPLE1" "$SAMPLE2"; do
  bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' \
    ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_filtered.vcf.gz \
    > ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_variants.bed
  
  bedtools genomecov \
    -ibam ${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}_recalibrated.bam \
    -bg \
    > ${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}_coverage.bedgraph
done
