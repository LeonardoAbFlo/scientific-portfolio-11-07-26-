#!/bin/bash

# Title: Variant Extraction and Counting
# Description: This script identifies and counts genetic variants (SNPs and Indels) 
#              from BAM files using BCFtools and BEDTools. It filters variants within
#              specific genes of interest defined in a BED file.
# Usage: Run this script from the terminal. Ensure that the necessary dependencies
#        (BCFtools, BEDTools, Samtools) are installed and accessible.
#        Example:
#        ./04_verify_variants_manually.sh
#
# Requirements:
# - BCFtools, BEDTools, and Samtools installed
# - Reference genome in FASTA format
# - Indexed BAM files with their respective BAI index
# - BED file specifying genes of interest

source $(conda info --base)/etc/profile.d/conda.sh

