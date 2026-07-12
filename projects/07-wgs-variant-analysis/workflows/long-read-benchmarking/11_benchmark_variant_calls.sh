#!/usr/bin/env bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate hap.py
hap.py \
    ${INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
    ${OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
    -f "${INPUT_DIR}/${BASELINE_BED_FILE_PATH}" \
    -r "${INPUT_DIR}/${REF}" \
    -o "${OUTPUT_DIR}/happy" \
    -l ${CONTIGS}:${START_POS}-${END_POS} \
    --engine=vcfeval \
    --threads="${THREADS}" \
    --pass-only