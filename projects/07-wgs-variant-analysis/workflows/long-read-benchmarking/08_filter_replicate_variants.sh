#!/bin/bash
# Filtering Variants Consistently Found in ≥2 Replicates

source $(conda info --base)/etc/profile.d/conda.sh

REFERENCE="/path/to/Desktop/benchmarking/references/tbdb.fasta"      # Path to your reference genome
INPUT_VCF="/path/to/Desktop/benchmarking/bcftools/test/combined_variants_v2.vcf.gz"
OUTPUT_VCF="/path/to/Desktop/benchmarking/bcftools/test/mtb_filtered_variants.vcf.gz"


# Define barcode groups per treatment
declare -A TREATMENTS=(
  ["70H"]="barcode08 barcode16 barcode24"
  ["80H"]="barcode07 barcode15 barcode23"
  ["85H"]="barcode06 barcode14 barcode22"
  ["90H"]="barcode05 barcode13 barcode21"
  ["95H"]="barcode04 barcode12 barcode20"
  ["97H"]="barcode03 barcode11 barcode19"
  ["99H"]="barcode02 barcode10 barcode18"
  ["100H"]="barcode01 barcode09 barcode17"
)

# Generate sample lists per treatment
for treatment in "${!treatments[@]}"; do
  echo "${treatments[$treatment]}" | tr ' ' '\n' > "${treatment}_samples.txt"
done

conda activate bcftools 

# Filter variants found in ≥2 out of 3 replicates within each treatment
bcftools view -i '
  (GT[@T1_samples.txt]="alt")>=2 ||
  (GT[@T2_samples.txt]="alt")>=2 ||
  (GT[@T3_samples.txt]="alt")>=2 ||
  (GT[@T4_samples.txt]="alt")>=2 ||
  (GT[@T5_samples.txt]="alt")>=2 ||
  (GT[@T6_samples.txt]="alt")>=2 ||
  (GT[@T7_samples.txt]="alt")>=2 ||
  (GT[@T8_samples.txt]="alt")>=2
' "${INPUT_VCF}" -Oz -o "${OUTPUT_VCF}"

tabix -p vcf "${OUTPUT_VCF}"

echo "Filtered VCF: ${OUTPUT_VCF}"