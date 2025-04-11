#!/bin/bash

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
# Relevant files and directories for analysis.

SRA=SRR14329362
REF_ID=NG_059281

RESULTS_DIR=BINF5002_project
RAW_DIR=${RESULTS_DIR}/rawdata
QC_REPORT=${RESULTS_DIR}/qc
TRIMMED_DIR=${RESULTS_DIR}/trimmed

echo Setting up directories...
mkdir -p $RESULTS_DIR
mkdir -p $RAW_DATA $QC_REPORT $TRIMMED_DIR
echo Completed setting up directories.

# --------------------------------------------------------------------------------------
# Download reference fasta sequence.
echo Downloading reference sequence...
efetch -db nucleotide -id $REF_ID -format fasta > "${RAW_DIR}/reference_${REF_ID}.fasta"
echo Completed downloading reference sequence.

# Download raw sequencing data.
echo Downloading FASTQ file...
prefetch $SRA -O $RAW_DIR 
fastq-dump ./${RAW_DIR}/${SRA} -O $RAW_DIR
echo Completed downloading raw sequencing data.

# Check if the raw FASTQ sequencing file has been downloaded.
if [ ! -s "$RAW_DIR/${SRA}.fastq" ]; then
    # If the user provided an invalid accession number, technically the efetch command still succeeded, but return an empty file.
    # -s test if a file existed and is bigger than 0 (empty file). So " ! -s <file> " returns TRUE when the file is empty or not exist. 
    echo "An error has occurred: The result is empty. Please check the accession number."
    exit
fi


# --------------------------------------------------------------------------------------
# We also need to check if the sequencing is single end or paired end.
# Paired end reads will contain the /1 /2 identifiers or something similar.
# However, our sequences in the FASTQ file did not have anything representing that.
# So it's likely that our reads are single end.
# This is important for choosing the right tools and flags for bioinformatics commands.

# First do QC check.
echo Performing QC check on raw FASTQ file...
fastqc "${RAW_DIR}/${SRA}.fastq" -o $QC_REPORT
echo QC check completed.



# Use fastp to filter low quality read. 
# For our purpose (genomics variant calling), the read quality should be much higher at ~30.
# Since this is single end read, we will use -i and -o flags (lowercase).
echo Trimming low quality reads...
fastp -i "${RAW_DIR}/${SRA}.fastq" -q 30 -o "$TRIMMED_DIR/trimmed_${SRA}.fastq"
echo Completed trimming low quality reads.


# echo Performing QC check on trimmed FASTQ reads...















# Use tree to check for correct files and directories.
tree $RESULTS_DIR
