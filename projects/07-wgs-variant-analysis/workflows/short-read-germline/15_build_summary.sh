source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

set -euo pipefail

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32

for SAMPLE in "$SAMPLE1" "$SAMPLE2"; do
    echo "========================================="
    echo "WGS Quality Control Summary for ${SAMPLE}"
    echo "========================================="
    echo

    echo "=== ALIGNMENT QUALITY ==="
    echo -n "Mapping Rate: "
    grep -A1 "FIRST_OF_PAIR" ${ROOT_W}/results/qc/${SAMPLE}/${SAMPLE}_alignment_summary.txt | tail -1 | cut -f7
    echo "  (Benchmark: >95%)"
    echo

    echo -n "Duplicate Rate: "
    grep -A1 "LIBRARY" ${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}_duplicate_metrics.txt | tail -1 | cut -f9
    echo "  (Benchmark: <30%)"
    echo

    echo -n "Mean Insert Size: "
    grep -A1 "MEDIAN_INSERT_SIZE" ${ROOT_W}/results/qc/${SAMPLE}/${SAMPLE}_insert_size_metrics.txt | tail -1 | cut -f1
    echo "bp  (Benchmark: 300-500bp)"
    echo

    echo "=== VARIANT CALLING QUALITY ==="
    echo -n "Total SNPs: "
    cat ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_snp_count.txt
    echo "  (Benchmark: 4-5 million)"
    echo

    echo -n "Total Indels: "
    cat ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_indel_count.txt
    echo "  (Benchmark: 0.5-0.8 million)"
    echo

    echo -n "Ti/Tv Ratio: "
    grep -v "^#" ${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_variant_stats.txt | grep "TSTV" | cut -f5
    echo "  (Benchmark: 2.0-2.1)"
    echo

    echo "For detailed variant annotation statistics, open:"
    echo "${ROOT_W}/results/var/${SAMPLE}/${SAMPLE}_annotation_stats.html"
    echo "========================================="
done
