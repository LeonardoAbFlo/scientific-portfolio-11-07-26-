source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project

cd $ROOT_W
mkdir -p matched
cd matched

SOMATIC_DIR=$ROOT_W/matched

mkdir -p {input_data,raw_calls,filtered_calls,contamination_analysis}
mkdir -p {converted_tables,maf_files,qc_reports}

cd input_data

TUMOR=SRR7890824
NORMAL=SRR7890827

for SAMPLE in "$TUMOR" "$NORMAL"; do
ln -s ${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}_recalibrated.bam ${SAMPLE}_recalibrated.bam
ln -s ${ROOT_W}/results/aligned/${SAMPLE}/${SAMPLE}_recalibrated.bai ${SAMPLE}_recalibrated.bai
done
