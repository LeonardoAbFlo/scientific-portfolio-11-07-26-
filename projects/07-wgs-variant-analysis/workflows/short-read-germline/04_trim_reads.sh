source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32

trim_galore \
  --paired \
  --quality 20 \
  --length 50 \
  --stringency 3 \
  --trim-n \
  --fastqc \
  --output_dir $ROOT_W/results/trimmed/${SAMPLE1} \
  $ROOT_W/data/${SAMPLE1}_1.fastq.gz \
  $ROOT_W/data/${SAMPLE1}_2.fastq.gz

trim_galore \
  --paired \
  --quality 20 \
  --length 50 \
  --stringency 3 \
  --trim-n \
  --fastqc \
  --output_dir $ROOT_W/results/trimmed/${SAMPLE2} \
  $ROOT_W/data/${SAMPLE2}_1.fastq.gz \
  $ROOT_W/data/${SAMPLE2}_2.fastq.gz
