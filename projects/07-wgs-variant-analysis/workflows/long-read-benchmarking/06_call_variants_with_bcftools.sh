#!/usr/bin/env bash

# Title: Variant Calling (Nanopore) with Pre-Sorted BAM
# Description:
#   This script calls variants using BCFtools (v1.21) on pre-sorted BAM files that were 
#   already generated from Oxford Nanopore Technologies (ONT) reads. If needed,
#   you can uncomment the alignment and sorting sections or place them behind a
#   conditional check to do everything in one pipeline.
#
# Usage:
#   1. Make the script executable:
#        chmod +x nanopore_variant_calling.sh
#   2. Run the script:
#        ./nanopore_variant_calling.sh
#   3. Ensure your sorted BAM files are in the directory specified by SORTED_BAM_DIR.
#
# Requirements:
#   - Conda or another environment with:
#       bcftools
#       samtools
#   - A reference genome in FASTA format (indexed with samtools faidx).
#   - Pre-sorted, indexed BAM files in a known directory.
#   - Sufficient disk space for output VCF files.

source $(conda info --base)/etc/profile.d/conda.sh

# ------------------------ USER-DEFINED VARIABLES ------------------------ #
REFERENCE="/path/to/Desktop/benchmarking/references/tbdb.fasta"      # Path to your reference genome
SORTED_BAM_DIR="/path/to/Desktop/benchmarking/depth_results/sorted_mapped_bam/standard" 
REGIONS="/path/to/Desktop/benchmarking/references/tbdb.bed"
OUTPUT_DIR="/path/to/Desktop/benchmarking/bcftools/test"       # Where final outputs (VCF, etc.) will be placed
THREADS=30
MAX_DEPTH=1500      # Adjust based on expected coverage
# ----------------------------------------------------------------------- #

mkdir -p "${OUTPUT_DIR}"

# (OPTIONAL) REFERENCE INDEXING
# conda activate samtools
# samtools faidx "${REFERENCE}"

# Define samples clearly (barcode-to-treatment mapping)
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

# Create a sample list file (for labeling in VCF)
SAMPLE_LIST="${OUTPUT_DIR}/sample_list.txt"
> "${SAMPLE_LIST}"  # empty or create file
for treatment in "${!TREATMENTS[@]}"; do
    for barcode in ${TREATMENTS[$treatment]}; do
        BAM_FILE="${SORTED_BAM_DIR}/${barcode}.sorted.bam"
        if [[ -f "${BAM_FILE}" ]]; then
            echo -e "${barcode}\t${treatment}" >> "${SAMPLE_LIST}"
        else
            echo "Warning: ${BAM_FILE} not found."
        fi
    done
done

echo "Sample\tTreatment"
cat "${SAMPLE_LIST}"

conda activate bcftools

# Joint variant calling with bacterial-specific parameters
ALL_BAMS=$(ls "${SORTED_BAM_DIR}"/barcode*.sorted.bam | tr '\n' ' ')

echo "Calling variants for all sorted BAMs in ${SORTED_BAM_DIR}..."

bcftools mpileup -f "${REFERENCE}" \
  -d ${MAX_DEPTH} \
  --threads 30 \
  --min-MQ 5 \
  --min-BQ 5 \
  -R ${REGIONS} \
  "${ALL_BAMS}" | 
  bcftools call -mv \
  --ploidy 1 \
  -Ov \
  --prior 0.1 \
  -P 20 \
  --keep-alts \
  -o "${OUTPUT_DIR}/test.vcf"

echo "Variant calling complete. Output at ${OUTPUT_DIR}/combined_variants_v2.vcf"

# OPTIONAL: Compress and index VCF for downstream analysis
bgzip -c "${OUTPUT_DIR}/combined_variants_v2.vcf" > "${OUTPUT_DIR}/combined_variants_v2.vcf.gz"
tabix -p vcf "${OUTPUT_DIR}/combined_variants_v2.vcf.gz"
