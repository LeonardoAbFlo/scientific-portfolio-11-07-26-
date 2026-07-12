#!/bin/bash

# Title: Depth
# Description: This script processes FASTQ files from ONT sequencing by mapping them to a reference genome,
#              converting and sorting BAM files, calculating depth, and extracting coverage information.
# Usage: Run this script from the terminal. Ensure that the necessary dependencies (Minimap2, Samtools, Conda)
#        are installed and configured properly.
#        Example:
#        ./depth.sh
#
# Requirements:
# - Conda environment with minimap2 and samtools installed
# - Reference genome in FASTA format
# - BED file defining target regions
# - FASTQ files stored in subdirectories within the FASTQ_DIR

source $(conda info --base)/etc/profile.d/conda.sh

BLUE="\033[1;34m"
GREEN="\033[1;32m"
NC="\033[0m"

# Define directory paths
PROJECT_DIR="/path/to/Desktop/benchmarking/first_run/depth_results"
REF_GENOME="/path/to/Desktop/benchmarking/references/tbdb.fasta"
BED_FILE="/path/to/Desktop/benchmarking/references/tbdb.bed"  # BED file for specific regions
FASTQ_DIR="/path/to/Desktop/benchmarking/first_run"
LOG_DIR="/path/to/Desktop/benchmarking/scripts/logs/first_run"
mkdir -p "$LOG_DIR" "$PROJECT_DIR"

SCRIPT_NAME=$(basename "$0")
LOG_NAME="${SCRIPT_NAME#[0-9]*-}"
LOG_NAME="${LOG_NAME%.sh}.log" 
LOG_PATH="${LOG_DIR}/${LOG_NAME}"

exec > >(tee -a "$LOG_PATH") 2>&1
start_time=$(date +%s)

echo -e "${BLUE}==================================================${NC}" | tee -a "$LOG_PATH"
echo -e "${BLUE}Script: $SCRIPT_NAME started at $(date)${NC}" | tee -a "$LOG_PATH"
echo -e "${BLUE}==================================================${NC}" | tee -a "$LOG_PATH"


# Process each subdirectory inside fastq/
for sample_dir in "$FASTQ_DIR"/*; do
    sample_name=$(basename "$sample_dir")
    echo "Processing sample folder: $sample_name..."

    # Define sample-specific directories
    MAPPED_DATA_DIR="$PROJECT_DIR/mapped_data/$sample_name"
    ORDERED_SAM_DIR="$PROJECT_DIR/ordered_sam/$sample_name"
    MAPPED_BAM_DIR="$PROJECT_DIR/mapped_bam/$sample_name"
    SORTED_MAPPED_BAM_DIR="$PROJECT_DIR/sorted_mapped_bam/$sample_name"
    DEPTH_DIR="$PROJECT_DIR/processed_data/depth/$sample_name"
    GENE_DEPTH_DIR="$PROJECT_DIR/processed_data/gene_depth/$sample_name"
    COVERAGE_DIR="$PROJECT_DIR/processed_data/coverage/$sample_name"

    # Create directories if they do not exist
    mkdir -p "$MAPPED_DATA_DIR" "$ORDERED_SAM_DIR" "$MAPPED_BAM_DIR" "$SORTED_MAPPED_BAM_DIR" "$DEPTH_DIR" "$GENE_DEPTH_DIR" "$COVERAGE_DIR"

    # Process each FASTQ file within the sample directory
    for fastq_file in "$sample_dir"/*.fastq; do
        file_name=$(basename "$fastq_file" .fastq)
        echo "-----------------------------------------------------------"
        echo "Processing file $file_name in sample folder $sample_name..."

        # Step 1: Map with minimap2
        conda activate minimap2
        sam_file="$MAPPED_DATA_DIR/$file_name.sam"
        minimap2 -t 26 -ax map-ont "$REF_GENOME" "$fastq_file" > "$sam_file"
        echo "Mapping completed for $file_name."

        # Step 2: Convert SAM to BAM
        conda activate samtools
        mapped_bam_file="$MAPPED_BAM_DIR/$file_name.unsorted.bam"
        samtools view -b "$sam_file" > "$mapped_bam_file"
        echo "BAM file generated for $file_name."

        # Step 3: Sort BAM file
        sorted_mapped_bam_file="$SORTED_MAPPED_BAM_DIR/$file_name.sorted.bam"
        samtools sort -@ 26 -o "$sorted_mapped_bam_file" "$mapped_bam_file"
        echo "Sorted BAM file generated for $file_name."
        samtools index "${sorted_mapped_bam_file}"
        

        # Step 4: Calculate genome-wide depth
        depth_file="$DEPTH_DIR/depth_$file_name.txt"
        samtools depth -@ 26 -a "$sorted_mapped_bam_file" > "$depth_file"
        echo "Genome-wide depth calculated for $file_name."

        # Step 5: Calculate depth for specific regions (from BED file)
        gene_depth_file="$GENE_DEPTH_DIR/gene_depth_$file_name.txt"
        samtools depth -@ 26 -b "$BED_FILE" "$sorted_mapped_bam_file" > "$gene_depth_file"
        echo "Region-specific depth calculated for $file_name."

        # Step 6: Calculate coverage
        coverage_file="$COVERAGE_DIR/coverage_$file_name.txt"
        samtools coverage "$sorted_mapped_bam_file" > "$coverage_file"
        echo "Coverage calculated for $file_name."

        echo -e "${GREEN}Processing completed for $file_name in sample folder $sample_name.${NC}"
    done

done
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

# Print completion message in blue
echo -e "${BLUE}==================================================${NC}" | tee -a "$LOG_PATH"
echo -e "${BLUE}Script: $SCRIPT_NAME completed at $(date)${NC}" | tee -a "$LOG_PATH"
echo -e "${BLUE}Total execution time: ${hours}h ${minutes}m ${seconds}s${NC}" | tee -a "$LOG_PATH"
echo -e "${BLUE}==================================================${NC}" | tee -a "$LOG_PATH"

# Save execution time separately
echo "Processing completed for all samples in fastq."



: '
# Example alignment/sorting process (COMMENTED OUT by default):

FASTQ_DIR="fastq_data"

find "${FASTQ_DIR}" -type f -name "*.fastq" | while read -r FASTQ_FILE; do
    SAMPLE_NAME=$(basename "${FASTQ_FILE}" .fastq)
    echo "Processing sample: ${SAMPLE_NAME}"

    # 2a. MAPPING with minimap2
    minimap2 -ax map-ont -t ${THREADS} "${REFERENCE}" "${FASTQ_FILE}" \
        | samtools view -bS - \
        > "${OUTPUT_DIR}/${SAMPLE_NAME}.unsorted.bam"

    # 2b. SORTING
    samtools sort -@ ${THREADS} -o "${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.bam" \
        "${OUTPUT_DIR}/${SAMPLE_NAME}.unsorted.bam"

    # 2c. INDEXING
    samtools index "${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.bam"

    rm "${OUTPUT_DIR}/${SAMPLE_NAME}.unsorted.bam"
done
'
