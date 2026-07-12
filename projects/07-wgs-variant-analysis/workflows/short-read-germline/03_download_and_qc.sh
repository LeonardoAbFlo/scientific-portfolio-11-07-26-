source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference/Homo_sapiens_assembly38.fasta"
SAMPLE1=SRR7890824
SAMPLE2=SRR7890827
THREADS=32

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR7890827/SRR7890827 --output-document=SRR7890827.sra  # normal sample
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR7890824/SRR7890824 --output-document=SRR7890824.sra  # tumor sample

fastq-dump --split-files $SAMPLE1.sra --gzip -O $ROOT_W/data/
fastq-dump --split-files $SAMPLE2.sra --gzip -O $ROOT_W/data/

mkdir -p $ROOT_W/results/{qc,trimmed,aligned,recal,var}/${SAMPLE1}
mkdir -p $ROOT_W/results/{qc,trimmed,aligned,recal,var}/${SAMPLE2}

fastqc \
    $ROOT_W/data/${SAMPLE1}_1.fastq.gz \
    $ROOT_W/data/${SAMPLE1}_2.fastq.gz \
    -o $ROOT_W/results/qc/${SAMPLE1}

fastqc \
    $ROOT_W/data/${SAMPLE2}_1.fastq.gz \
    $ROOT_W/data/${SAMPLE2}_2.fastq.gz \
    -o $ROOT_W/results/qc/${SAMPLE2}
