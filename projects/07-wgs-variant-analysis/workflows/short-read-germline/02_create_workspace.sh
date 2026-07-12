source ~/anaconda3/etc/profile.d/conda.sh
conda activate wgs

ROOT_W=/path/to/wgs_project

mkdir -p $ROOT_W/{reference,data,results}
cd $ROOT_W/reference

wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
bwa index Homo_sapiens_assembly38.fasta

wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi
