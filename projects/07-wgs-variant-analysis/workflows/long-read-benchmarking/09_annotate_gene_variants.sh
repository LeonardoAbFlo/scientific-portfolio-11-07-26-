#!/bin/bash

FILTERED_VCF="mtb_filtered_variants.vcf.gz"
GENE_BED="resistance_genes.bed"
ANNOTATED_OUTPUT="mtb_resistance_variants.tsv"

# Intersect variants with resistance genes
bedtools intersect \
  -a <(bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\t[%SAMPLE:%GT\t]\n' "${FILTERED_VCF}") \
  -b "${GENE_BED}" \
  -wa -wb > "${ANNOTATED_OUTPUT}"

echo "Annotated variants saved to: ${ANNOTATED_OUTPUT}"