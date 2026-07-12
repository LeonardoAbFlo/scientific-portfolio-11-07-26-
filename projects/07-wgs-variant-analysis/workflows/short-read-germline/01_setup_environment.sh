source ~/anaconda3/etc/profile.d/conda.sh
conda create -y -n wgs python=3.9
conda activate wgs
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install -y gatk4 bwa samtools picard trim-galore fastqc bcftools snpeff bedtools vcftools tabix tree sra-toolkit ncurses
conda list | egrep "^(gatk4|bwa|samtools|picard|trim-galore|fastqc|bcftools|snpeff|bedtools|vcftools|tabix|tree)\b"
