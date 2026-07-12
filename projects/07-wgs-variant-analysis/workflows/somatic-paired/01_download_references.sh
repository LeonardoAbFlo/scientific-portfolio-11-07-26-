ROOT_W=/path/to/wgs_project
REFERENCE="$ROOT_W/reference"

cd $REFERENCE
mkdir somatic_resources
cd somatic_resources

wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi

ln -s $REFERENCE/Homo_sapiens_assembly38.fasta ./
ln -s $REFERENCE/Homo_sapiens_assembly38.fasta.fai ./
ln -s $REFERENCE/Homo_sapiens_assembly38.dict ./
