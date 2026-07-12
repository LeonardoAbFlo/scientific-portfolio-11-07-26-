source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project
SOMATIC_DIR="${ROOT_W}/matched"

REFERENCE="${ROOT_W}/reference/somatic_resources/Homo_sapiens_assembly38.fasta"
GERMLINE_RESOURCE="${ROOT_W}/reference/somatic_resources/af-only-gnomad.hg38.vcf.gz"
PON="${ROOT_W}/reference/somatic_resources/1000g_pon.hg38.vcf.gz"

INPUT_DIR="${SOMATIC_DIR}/input_data"
OUTPUT_DIR="${SOMATIC_DIR}/raw_calls"

TUMOR=SRR7890824
NORMAL=SRR7890827

mkdir -p "$OUTPUT_DIR"

gatk Mutect2 \
    -R "$REFERENCE" \
    -I "${INPUT_DIR}/${TUMOR}_recalibrated.bam" \
    -I "${INPUT_DIR}/${NORMAL}_recalibrated.bam" \
    -tumor "$TUMOR" \
    -normal "$NORMAL" \
    --germline-resource "$GERMLINE_RESOURCE" \
    --panel-of-normals "$PON" \
    --f1r2-tar-gz "${OUTPUT_DIR}/${TUMOR}_f1r2.tar.gz" \
    -O "${OUTPUT_DIR}/${TUMOR}_raw.vcf.gz" \
    --native-pair-hmm-threads 32 \
    --max-reads-per-alignment-start 50

bcftools stats "${OUTPUT_DIR}/${TUMOR}_raw.vcf.gz" > "${OUTPUT_DIR}/${TUMOR}_raw_stats.txt"
