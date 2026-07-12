source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project
SOMATIC_DIR="${ROOT_W}/matched"

REFERENCE="${ROOT_W}/reference/somatic_resources/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SOMATIC_DIR}/raw_calls"
CONTAM_DIR="${SOMATIC_DIR}/contamination_analysis"

TUMOR=SRR7890824
NORMAL=SRR7890827

mkdir -p "${SOMATIC_DIR}/filtered_calls"

gatk LearnReadOrientationModel \
    -I "${OUTPUT_DIR}/${TUMOR}_f1r2.tar.gz" \
    -O "${OUTPUT_DIR}/${TUMOR}_orientation_model.tar.gz"

gatk FilterMutectCalls \
    -R "$REFERENCE" \
    -V "${OUTPUT_DIR}/${TUMOR}_raw.vcf.gz" \
    --contamination-table "${CONTAM_DIR}/${TUMOR}_contamination.table" \
    --tumor-segmentation "${CONTAM_DIR}/${TUMOR}_segments.table" \
    --ob-priors "${OUTPUT_DIR}/${TUMOR}_orientation_model.tar.gz" \
    -O "${SOMATIC_DIR}/filtered_calls/${TUMOR}_filtered.vcf.gz"

bcftools view -H "${SOMATIC_DIR}/filtered_calls/${TUMOR}_filtered.vcf.gz" \
    | cut -f7 \
    | sort \
    | uniq -c \
    | sort -nr
