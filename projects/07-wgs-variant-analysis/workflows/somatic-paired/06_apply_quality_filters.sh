source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project
SOMATIC_DIR="${ROOT_W}/matched"
OUTPUT_DIR="${SOMATIC_DIR}/raw_calls"

TUMOR=SRR7890824
NORMAL=SRR7890827

bcftools view -f PASS \
    "${SOMATIC_DIR}/filtered_calls/${TUMOR}_filtered.vcf.gz" \
    -O z \
    -o "${SOMATIC_DIR}/filtered_calls/${TUMOR}_pass.vcf.gz"

bcftools index -t "${SOMATIC_DIR}/filtered_calls/${TUMOR}_pass.vcf.gz"

bcftools filter \
    -i 'FORMAT/AF[0:0] >= 0.05 && FORMAT/DP[0:0] >= 10 && INFO/TLOD >= 6.3 && (FORMAT/AF[0:1] <= 0.03 || FORMAT/AF[0:1] == ".")' \
    "${SOMATIC_DIR}/filtered_calls/${TUMOR}_pass.vcf.gz" \
    -O z \
    -o "${SOMATIC_DIR}/filtered_calls/${TUMOR}_high_confidence.vcf.gz"

bcftools index -t "${SOMATIC_DIR}/filtered_calls/${TUMOR}_high_confidence.vcf.gz"

raw_count=$(bcftools view -H "${OUTPUT_DIR}/${TUMOR}_raw.vcf.gz" | wc -l)
pass_count=$(bcftools view -H "${SOMATIC_DIR}/filtered_calls/${TUMOR}_pass.vcf.gz" | wc -l)
hc_count=$(bcftools view -H "${SOMATIC_DIR}/filtered_calls/${TUMOR}_high_confidence.vcf.gz" | wc -l)

echo ""
echo "Final filtering cascade results:"
echo "  Raw Mutect2 calls: $raw_count"
echo "  PASS calls: $pass_count"
echo "  High-confidence calls: $hc_count"

if [ "$raw_count" -gt 0 ]; then
    echo "  Final success rate: $(echo "scale=2; $hc_count * 100 / $raw_count" | bc)%"
else
    echo "  Final success rate: NA (raw_count=0)"
fi
