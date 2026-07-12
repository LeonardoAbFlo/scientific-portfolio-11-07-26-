#!/bin/bash

# Title: Variant Verification
# Description: This script verifies the number of variants in VCF files generated from ONT sequencing data.
#              It processes each sample directory, extracts variant statistics, and logs the results.
# Usage: Run this script from the terminal. Ensure that bcftools and conda are properly installed.
#        Example:
#        ./variant_verification.sh
#
# Requirements:
# - Conda environment with bcftools installed
# - VCF files for variant analysis
# - Reference genome in FASTA format
# - Indexed VCF files (.vcf.gz and .tbi)

source $(conda info --base)/etc/profile.d/conda.sh

BLUE="\033[1;34m"
GREEN="\033[1;32m"
NC="\033[0m"

# Define directory paths
PROJECT_DIR="/path/to/Desktop/benchmarking/variant_results"
VCF_DIR="/path/to/Desktop/benchmarking/vcf"
LOG_DIR="/path/to/Desktop/benchmarking/scripts/logs"
SCRIPT_NAME=$(basename "$0")
LOG_NAME="${SCRIPT_NAME#[0-9]*-}"
LOG_NAME="${LOG_NAME%.sh}.log"
LOG_PATH="${LOG_DIR}/${LOG_NAME}"

exec > >(tee -a "$LOG_PATH") 2>&1
start_time=$(date +%s)

echo -e "${BLUE}==================================================${NC}" | tee -a "$LOG_PATH"
echo -e "${BLUE}Script: $SCRIPT_NAME started at $(date)${NC}" | tee -a "$LOG_PATH"
echo -e "${BLUE}==================================================${NC}" | tee -a "$LOG_PATH"

# Create output directory if it does not exist
mkdir -p "$PROJECT_DIR"

# Process each VCF file
for vcf_file in "$VCF_DIR"/*.vcf.gz; do
    sample_name=$(basename "$vcf_file" .vcf.gz)
    echo "Processing VCF file: $sample_name..."

    # Define output file
    variant_stats_file="$PROJECT_DIR/${sample_name}_variant_stats.txt"

    # Step 1: Check if VCF file is indexed
    if [ ! -f "${vcf_file}.tbi" ]; then
        echo "Indexing VCF file: $sample_name..."
        conda activate bcftools
        bcftools index "$vcf_file"
    fi

    # Step 2: Count total variants
    total_variants=$(bcftools view -H "$vcf_file" | wc -l)
    echo "Total variants in $sample_name: $total_variants"

    # Step 3: Count SNPs and INDELs separately
    snps=$(bcftools view -v snps "$vcf_file" | wc -l)
    indels=$(bcftools view -v indels "$vcf_file" | wc -l)
    echo "SNPs in $sample_name: $snps"
    echo "INDELs in $sample_name: $indels"

    # Step 4: Count high-quality variants (QUAL ≥ 30)
    high_qual_variants=$(bcftools view -i 'QUAL>=30' "$vcf_file" | wc -l)
    echo "High-quality variants (QUAL >= 30) in $sample_name: $high_qual_variants"

    # Step 5: Extract basic statistics using bcftools stats
    bcftools stats "$vcf_file" > "$variant_stats_file"
    echo "Variant statistics saved to $variant_stats_file"

    echo -e "${GREEN}Processing completed for $sample_name.${NC}"
done

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo -e "${BLUE}==================================================${NC}" | tee -a "$LOG_PATH"
echo -e "${BLUE}Script: $SCRIPT_NAME completed at $(date)${NC}" | tee -a "$LOG_PATH"
echo -e "${BLUE}Total execution time: ${hours}h ${minutes}m ${seconds}s${NC}" | tee -a "$LOG_PATH"
echo -e "${BLUE}==================================================${NC}" | tee -a "$LOG_PATH"

echo "Variant verification completed for all samples."