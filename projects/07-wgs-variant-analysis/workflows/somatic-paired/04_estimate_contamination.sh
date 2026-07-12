source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project
SOMATIC_DIR="${ROOT_W}/matched"

INPUT_DIR="${SOMATIC_DIR}/input_data"
CONTAM_DIR="${SOMATIC_DIR}/contamination_analysis"
COMMON_VARIANTS="${ROOT_W}/reference/somatic_resources/small_exac_common_3.hg38.vcf.gz"

TUMOR=SRR7890824
NORMAL=SRR7890827

mkdir -p "$CONTAM_DIR"

gatk GetPileupSummaries \
    -I "${INPUT_DIR}/${TUMOR}_recalibrated.bam" \
    -V "$COMMON_VARIANTS" \
    -L "$COMMON_VARIANTS" \
    -O "${CONTAM_DIR}/${TUMOR}_pileups.table"

gatk GetPileupSummaries \
    -I "${INPUT_DIR}/${NORMAL}_recalibrated.bam" \
    -V "$COMMON_VARIANTS" \
    -L "$COMMON_VARIANTS" \
    -O "${CONTAM_DIR}/${NORMAL}_pileups.table"

gatk CalculateContamination \
    -I "${CONTAM_DIR}/${TUMOR}_pileups.table" \
    -matched "${CONTAM_DIR}/${NORMAL}_pileups.table" \
    -O "${CONTAM_DIR}/${TUMOR}_contamination.table" \
    --tumor-segmentation "${CONTAM_DIR}/${TUMOR}_segments.table"
