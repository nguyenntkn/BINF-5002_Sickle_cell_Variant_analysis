
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
}

# usage 
