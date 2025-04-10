# 

function usage {
    echo Project: BINF5002 final project - Genomics variant analysis of sickle cell anemia
    echo Authors: Nikita Ch, Simran, Pragati Dwivedi, Tran Khanh Nguyen Nguyen
    echo Description:
    echo Bioconda packages required: entrez-direct, sra-tools fastqc, fastp, cutadapt, trimmomatic, multiqc, samtools...
    echo 
    echo NOTE: This script is designed for working with restricted licensing within Cocalc environment. Some manual set up is required before running the script.
    echo Set up a conda environment:
    echo anaconda2023
    echo conda create --prefix ./bioenv
    echo conda activate /home/user/bioenv
    echo conda install -n bioconda entrez-direct sra-tools fastqc fastp cutadapt trimmomatic multiqc samtools -y
    exit 1
}

# usage 

# --------------------------------------------------------------------------------------
# Relevant file

SRA_RAW=SRR14329362
REF_ID=NG_059281

RESULTS_DIR=BINF5002_project
RAW_DATA=${RESULTS_DIR}/rawdata


echo Setting up directories...
mkdir -p $RESULTS_DIR
mkdir -p $RAW_DATA
echo Completed setting up directories.

# --------------------------------------------------------------------------------------
# Download reference fasta sequence
echo Downloading reference sequence...
efetch -db nucleotide -id $REF_ID -format fasta > ${RAW_DATA}/reference_${REF_ID}.fasta
echo Completed downloading reference sequence.

# Download raw sequencing data
echo Downloading FASTQ file...
prefetch $SRA_RAW -O $RAW_DATA 
fastq-dump ./${RAW_DATA}/${SRA_RAW} -O $RAW_DATA
echo Completed downloading raw sequencing data.

tree $RESULTS_DIR


