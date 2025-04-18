#!/bin/bash

function usage {
    echo Project: BINF5002 final project - Genomics variant analysis of sickle cell anemia
    echo Authors: Nikita Chaudhari, Simran Wilasra, Pragati Dwivedi, Tran Khanh Nguyen Nguyen
    echo Description:
    echo Bioconda packages required: entrez-direct, sra-tools, fastqc, fastp, cutadapt, trimmomatic, multiqc, samtools...
    echo Actually used: entrez-direct, sra-tools, fastp, bwa
    echo 
    echo NOTE: This script is designed for working with restricted licensing within Cocalc environment. Some manual set up is required before running the script.
    echo Set up a conda environment:
    echo anaconda2023
    echo conda create --prefix ./bioenv
    echo conda activate /home/user/bioenv
    echo conda install -c bioconda entrez-direct sra-tools fastqc fastp cutadapt trimmomatic multiqc samtools -y
    echo conda install bwa
    echo conda install -c bioconda gatk4 snpEff -y
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
ALIGNED_DIR=${RESULTS_DIR}/aligned
VARIANT_DIR=${RESULTS_DIR}/variants
SNPEFF_DIR=${RESULTS_DIR}/snpEff
SNPEFF_DATA_DIR=${SNPEFF_DIR}/data/reference_dbd

echo Setting up directories...
mkdir -p $RESULTS_DIR
mkdir -p $RAW_DIR $QC_REPORT $TRIMMED_DIR $ALIGNED_DIR
mkdir -p $VARIANT_DIR
mkdir -p $SNPEFF_DATA_DIR
mkdir -p $ANNOTATED_DIR
mkdir -p $SNPEFF_DIR $SNPEFF_DATA_DIR
echo Completed setting up directories.

# --------------------------------------------------------------------------------------
# DOWNLOADING REFERENCE SEQUENCE AND RAW SEQUENCING SEQUENCES

# Download reference fasta sequence.
echo Downloading reference sequence...
efetch -db nucleotide -id $REF_ID -format fasta > "${RAW_DIR}/reference_${REF_ID}.fasta"
echo Completed downloading reference sequence.

# Check if the reference sequence file has been downloaded.
if [ ! -s "${RAW_DIR}/reference_${REF_ID}.fasta" ]; then
    # If the user provided an invalid accession number, technically the efetch command still succeeded, but return an empty file.
    # -s test if a file existed and is bigger than 0 (empty file). So " ! -s <file> " returns TRUE when the file is empty or not exist.     
    echo An error has occurred: The reference FASTA file is empty. Please check the accession number.
    exit 1
fi

# Download raw sequencing data.
echo Downloading FASTQ file...
prefetch $SRA -O $RAW_DIR 
fastq-dump ./${RAW_DIR}/${SRA} -O $RAW_DIR
echo Completed downloading raw sequencing data.

# Check if the raw FASTQ sequencing file has been downloaded.
if [ ! -s "${RAW_DIR}/${SRA}.fastq" ]; then
    echo An error has occurred: The FASTQ file is empty. Please check the accession number.
    exit 1
fi


# --------------------------------------------------------------------------------------
# QC CHECK 

# We also need to check if the sequencing is single end or paired end.
# Paired end reads will contain the /1 /2 identifiers or something similar.
# However, our sequences in the FASTQ file did not have anything representing that.
# So it's likely that our reads are single end.
# This is important for choosing the right tools and flags for bioinformatics commands.


# # Instead of creating a pipeline with fastqc and cutadapt/ trimmomatic, we can just use fastp
# # First do QC check.
# echo Performing QC check on raw FASTQ file...
# fastqc "${RAW_DIR}/${SRA}.fastq" -o $QC_REPORT
# echo QC check completed.


# Use fastp to filter low quality read. 
# For our purpose (genomics variant calling), the read quality should be much higher at ~30.
# Since this is single end read, we will use -i and -o flags (lowercase).
echo Trimming low quality reads...
fastp -i "${RAW_DIR}/${SRA}.fastq" -o "${TRIMMED_DIR}/trimmed_${SRA}.fastq" -h "${QC_REPORT}/${SRA}_fastp_report.html" -q 30
echo Completed trimming low quality reads.


# --------------------------------------------------------------------------------------
# INDEXING REFERENCE SEQUENCE

# By default, bwa index creates several index files in the same directory as the input FASTA file, and it uses the same base name.
echo Indexing reference sequence...
bwa index "${RAW_DIR}/reference_${REF_ID}.fasta"
echo Completed indexing reference sequence.

# Sequence alignment
echo Aligning sequences...
bwa mem "${RAW_DIR}/reference_${REF_ID}.fasta" "${TRIMMED_DIR}/trimmed_${SRA}.fastq" > "${ALIGNED_DIR}/aligned_${SRA}.sam"
echo Completed aligning sequences.

#bwa mem with read group
#bwa mem -R "@RG\tID:SRR14329362\tLB:lib1\tPL:illumina\tPU:SRR14329362\tSM:sample1" "${RAW_DIR}/reference_${REF_ID}.fasta" "${TRIMMED_DIR}/trimmed_${SRA}.fastq" > "${ALIGNED_DIR}/aligned_${SRA}.sam"
#echo Inserting read groups  

#Convert SAM to sorted BAM; from human readable to binary file
echo Converting .SAM to sorted BAM file...
samtools view -b ${ALIGNED_DIR}/aligned_${SRA}.sam | samtools sort -o ${ALIGNED_DIR}/aligned_${SRA}.bam

#Install GATK - conda install -c bioconda gatk4 -y
#Validate BAM file
echo Validating BAM file..
gatk ValidateSamFile -I ${ALIGNED_DIR}/aligned_${SRA}.bam -MODE SUMMARY
#Error: missing read group -> therefore need to do bwa mem -R
#Now no errors found!

#Convert FASTQ to FASTA for BLAST
echo Converting .fastq to .fasta
seqtk seq -A ${TRIMMED_DIR}/trimmed_${SRA}.fastq > ${TRIMMED_DIR}/trimmed_${SRA}.fasta

#BLAST; ; can edit evalue if we're not getting the results we want; but we are 
echo Performing BLAST...
makeblastdb -in ${RAW_DIR}/reference_${REF_ID}.fasta -dbtype nucl -out blast/reference_${REF_ID}_db
blastn -query ${TRIMMED_DIR}/trimmed_${SRA}.fasta -db blast/reference_${REF_ID}_db -out blast/blast_output.txt -outfmt 0 -evalue 1e-10


#Indexing with samtools
echo Indexing reference genome with samtools...
samtools faidx $RAW_DIR/reference_${REF_ID}.fasta

#Creating FASTA dictionary using GATK
echo Creating FASTA dictionary using GATK
gatk CreateSequenceDictionary -R $RAW_DIR/reference_${REF_ID}.fasta -O $RAW_DIR/reference_${REF_ID}.dict

#Mark duplicates...
echo Marking nucleotide duplicates...
gatk MarkDuplicates -I ${ALIGNED_DIR}/aligned_${SRA}.bam -O ${ALIGNED_DIR}/deduplicated_${SRA}.bam -M ${ALIGNED_DIR}/${SRA}_duplicate_metrics.txt

#Index the deduplicated BAM file; creates .bai file to allow rapid access to regions of the BAM file
echo Indexing deduduplicated BAM file...
samtools index ${ALIGNED_DIR}/deduplicated_${SRA}.bam

#Calling variants... HaplotypeCaller detects SNPs and indels from the BAM file; outputs to .vcf file (text file)
echo Calling variants...
gatk HaplotypeCaller -R ${RAW_DIR}/reference_${REF_ID}.fasta -I ${ALIGNED_DIR}/deduplicated_${SRA}.bam -O ${VARIANT_DIR}/raw_variants.vcf

#Filter Variants
echo Filtering variants...
gatk VariantFiltration -R $RAW_DIR/reference_${REF_ID}.fasta -V $VARIANT_DIR/raw_variants.vcf -O $VARIANT_DIR/filtered_variants.vcf --filter-expression "QD < 2.0 || FS > 60.0" --filter-name FILTER
filtered variants same as raw -> because our disease is an SNP


# echo Downloading reference GenBank file for snpEff...
# efetch -db nucleotide -id $REF_ID -format genbank > $SNPEFF_DATA_DIR/HBB.gbk
# echo Downloaded GenBank file for snpEff!


# #Create snpEff condiguration file; clean up later by making new directory for config files 
# echo Creating custom snpEff configuration file...
# cat <<EOF > $SNPEFF_DIR/snpEff.config
# # Custom snpEff config for reference_db
# reference_db.genome : reference_db
# reference_db.fa : $(readlink -f $RAW_DIR/reference_${REF_ID}.fasta)
# reference_db.genbank : $(readlink -f $SNPEFF_DATA_DIR/HBB.gbk)
# EOF

# #conda install -c bioconda snpEff -y
# echo Building snpEff database...
# snpEff build -c $SNPEFF_DIR/snpEff.config -genbank -v -noCheckProtein reference_db
# echo Built snpEff database!

# #Exporting snpEff database (this is optional, we're just dumping the contents of the snpEff file into a text file for debugging)
# echo Exporting snpEff database...
# snpEff dump -c $SNPEFF_DIR/snpEff.config reference_db > $SNPEFF_DIR/snpEff_reference_db.txt 
# echo Exported snpEff database!


# echo Annotating variants with snpEff...
# snpEff -c $SNPEFF_DIR/snpEff.config -stats $SNPEFF_DIR/snpEff.html reference_db $VARIANT_DIR/filtered_variants.vcf > $ANNOTATED_DIR/annotated_variants.vcf
# echo Annotated variants with snpEff!


# Use tree to check for correct files and directories.
tree $RESULTS_DIR
