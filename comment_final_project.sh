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
    echo conda install -c bioconda entrez-direct sra-tools fastqc fastp multiqc samtools gatk4 snpEff -y
    echo conda install bwa seqtk -y
    echo Additional: conda install -c bioconda gatk4 snpEff cutadapt trimmomatic -y
    exit 1
}

# usage 

# ----------------------------PART 1: GETTING SET UP----------------------------------------------------------
# 
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
SNPEFF_DATA_DIR=${SNPEFF_DIR}/data/reference_db
BLAST=${RESULTS_DIR}/blast
ANNOTATED_DIR=${RESULTS_DIR}/annotated

echo SETTING UP DIRECTORIES...
mkdir -p $RESULTS_DIR
mkdir -p $RAW_DIR $QC_REPORT $TRIMMED_DIR $ALIGNED_DIR $VARIANT_DIR $ANNOTATED_DIR $SNPEFF_DIR $SNPEFF_DATA_DIR
mkdir -p $BLAST
# echo Completed setting up directories.

# Log output
exec > >(tee -i $RESULTS_DIR/log.txt)
exec 2>&1

# -------------------------PART 2: DOWNLOADING REFERENCE SEQUENCE AND RAW SEQUENCING SEQUENCES-----------------------------------

# Download reference fasta sequence.
echo DOWNLOADING REFERENCE SEQUENCE...
efetch -db nucleotide -id $REF_ID -format fasta > "${RAW_DIR}/reference_${REF_ID}.fasta"
# echo Completed downloading reference sequence.

# Check if the reference sequence file has been downloaded.
if [ ! -s "${RAW_DIR}/reference_${REF_ID}.fasta" ]; then
    # If the user provided an invalid accession number, technically the efetch command still succeeded, but return an empty file.
    # -s test if a file existed and is bigger than 0 (empty file). So " ! -s <file> " returns TRUE when the file is empty or not exist.     
    echo An error has occurred: The reference FASTA file is empty. Please check the accession number.
    exit 1
fi

# Download raw sequencing data.
echo DOWNLOADING RAW FASTQ FILE...
prefetch $SRA -O $RAW_DIR 
fastq-dump ./${RAW_DIR}/${SRA} -O $RAW_DIR
# echo Completed downloading raw sequencing data.

# Check if the raw FASTQ sequencing file has been downloaded.
if [ ! -s "${RAW_DIR}/${SRA}.fastq" ]; then
    echo An error has occurred: The FASTQ file is empty. Please check the accession number.
    exit 1
fi


# -----------------------------PART 3: QUALITY CONTROL AND DATA CLEANING-------------------------------------------


# We also need to check if the sequencing is single end or paired end.
# Paired end reads will contain the /1 /2 identifiers or something similar.
# However, our sequences in the FASTQ file did not have anything representing that.
# So it's likely that our reads are single end.
# This is important for choosing the right tools and flags for bioinformatics commands.


# Instead of creating a pipeline with fastqc and cutadapt/ trimmomatic, we can just use fastp
# Use fastp to filter low quality read. 
# For our purpose (genomics variant calling), the read quality should be much higher at ~30.
# Since this is single end read, we will use -i and -o flags (lowercase).
echo TRIMMING LOW QUALITY READS...
fastp -i "${RAW_DIR}/${SRA}.fastq" -o "${TRIMMED_DIR}/trimmed_${SRA}.fastq" -h "${QC_REPORT}/${SRA}_fastp_report.html" -q 30
# echo Completed trimming low quality reads.


# -----------------------PART 4: SEQUENCE ALIGNMENT AND VALIDATION----------------------------------------------------

# By default, bwa index creates several index files in the same directory as the input FASTA file, and it uses the same base name.
# NOTE: this is essential for bwa mem sequence alignment, and is different from samtools faidx.
echo INDEXING REFERENCE SEQUENCE FOR ...
bwa index "${RAW_DIR}/reference_${REF_ID}.fasta"
echo Completed indexing reference sequence.

# Sequence alignment
echo Aligning sequences...
bwa mem "${RAW_DIR}/reference_${REF_ID}.fasta" "${TRIMMED_DIR}/trimmed_${SRA}.fastq" > "${ALIGNED_DIR}/aligned_${SRA}.sam"
echo Completed aligning sequences.

# Adding read group to SAM file using bwa mem
echo INSERTING READ GROUP TO SAM FILE...  
bwa mem -R "@RG\tID:${SRA}\tLB:lib1\tPL:illumina\tPU:${SRA}\tSM:sample1" "${RAW_DIR}/reference_${REF_ID}.fasta" "${TRIMMED_DIR}/trimmed_${SRA}.fastq" > "${ALIGNED_DIR}/aligned_${SRA}.sam"

# bwa mem -R "@RG\tID:SRR14329362\tLB:Fig2e_EDFig8e_untreated_BoneMarrow_Donor1_mouse2\tPL:illumina\tPU:SRR14329362\tSM:Fig2e_EDFig8e_untreated_BoneMarrow_Donor1_mouse2" \
#     "${RAW_DIR}/reference_${REF_ID}.fasta" \
#     "${TRIMMED_DIR}/trimmed_${SRA}.fastq" > "${ALIGNED_DIR}/aligned_${SRA}.sam"
# @RG info can be found on SRA database with the SRA accession number
# ID: This is usually a unique identifier for the sequencing run. In this case, it could be the SRR number (e.g., SRR14329362), which is linked to your specific run.
# Example: SRR14329362
# LB (Library): This should refer to the name or identifier of the library that was prepared for sequencing. Based on the metadata you provided, the library name is:
# Example: Fig2e_EDFig8e_untreated_BoneMarrow_Donor1_mouse2
# PL (Platform): This refers to the sequencing platform used for sequencing the samples. Based on the information you provided, the sequencing platform is Illumina MiSeq.
# Example: illumina
# PU (Platform Unit): This typically refers to a unique identifier of the sequencing run, often the same as the sequencing run identifier (in this case, SRR14329362).
# Example: SRR14329362
# SM (Sample): This is the biological sample being sequenced. From the sample name provided, it appears to be:
# Example: Fig2e_EDFig8e_untreated_BoneMarrow_Donor1_mouse2

#Convert SAM to sorted BAM (from human readable to binary file)
echo CONVERTING SAM FILE TO SORTED BAM FILE...
samtools view -b ${ALIGNED_DIR}/aligned_${SRA}.sam | samtools sort -o ${ALIGNED_DIR}/aligned_${SRA}.bam

#Validate BAM file
echo VALIDATING BAM FILE...
gatk ValidateSamFile -I ${ALIGNED_DIR}/aligned_${SRA}.bam -MODE SUMMARY
#Error: missing read group -> therefore need to do bwa mem -R
#Now no errors found!
# Warning:	ValidateSamFile	NM validation cannot be performed without the reference. All other validations will still occur.
# NM is the edit distance (number of differences) from the reference sequence. We don't have to do this rn.

# ---------------------PART 5: BLAST (Why?)-------------------------------------------------------------

# Convert FASTQ to FASTA for BLAST
echo CONVERTING TRIMMED FASTQ FILE TO FASTA FORMAT FOR BLAST...
seqtk seq -A ${TRIMMED_DIR}/trimmed_${SRA}.fastq > ${TRIMMED_DIR}/trimmed_${SRA}.fasta

# #BLAST; ; can edit evalue if we're not getting the results we want; but we are 
echo PERFORMING BLAST...
makeblastdb -in ${RAW_DIR}/reference_${REF_ID}.fasta -dbtype nucl -out ${BLAST}/reference_${REF_ID}_db
blastn -query ${TRIMMED_DIR}/trimmed_${SRA}.fasta -db ${BLAST}/reference_${REF_ID}_db -out ${BLAST}/blast_output_${SRA}.txt -outfmt 0 -evalue 1e-10
# This will perform alignment on ALL sequences within the FASTQ file (converted to FASTA) so the blast_output_${SRA}.txt file will contain all of those alignments
# Blast returns many indexing files:
#   .ndb: Primary database file
#   .nhr: Header file    -> faster lookups
#   .nin: Index file     -> faster lookups
#   .nsq: Sequence file  -> faster lookups
#   .ntk: FAT index file (filtering and translation)
#   .ntf: Optional for translation offset information (when translated into different reading frames)
#   .not: optional for nulti-threaded database
# But the actual alignment is stored in ${BLAST}/blast_output_${SRA}.txt file

# -----------------------------------PART 6: POST ALIGNMENT PROCESSING---------------------------------------------

# Indexing with samtools
# Samtools faidx is different from bwa index. Bwa index creates files important for sequence alignment,
# while samtools faidx simply creates a .fai file, which is relevant for downstream analysis such as gatk tools.
echo INDEXING REFERENCE SEQUENCE FOR DOWNSTREAM VARIANT CALLING ANALYSIS...
samtools faidx $RAW_DIR/reference_${REF_ID}.fasta

# Creating FASTA dictionary using GATK
# GATK needs to know the order and names of the contigs (like chromosomes) in the reference sequence, which it reads from the .dict file. 
# This is essential for things like sorting and variant calling.
echo CREATING GATK DICTIONARY FROM REFERENCE FASTA FILE
gatk CreateSequenceDictionary -R $RAW_DIR/reference_${REF_ID}.fasta -O $RAW_DIR/reference_${REF_ID}.dict

# Mark duplicates.
# Identifies duplicate reads from PCR amplifications in sorted BAM files. 
# Flag duplicates to prevent it from biasing variant calling, optional removal of duplicates.
# Generates metrics to assess duplication levels.
echo MARKING NUCLEOTIDE DUPLICATES...
gatk MarkDuplicates -I ${ALIGNED_DIR}/aligned_${SRA}.bam -O ${ALIGNED_DIR}/deduplicated_${SRA}.bam -M ${ALIGNED_DIR}/${SRA}_duplicate_metrics.txt

# Index the deduplicated BAM file; creates .bai file to allow rapid access to regions of the BAM file
echo INDEXING DEDUPLICATED FILE...
samtools index ${ALIGNED_DIR}/deduplicated_${SRA}.bam

# Calling variants... HaplotypeCaller detects SNPs and indels from the BAM file; outputs to .vcf file (text file)
echo CALLING VARIANTS...
gatk HaplotypeCaller -R ${RAW_DIR}/reference_${REF_ID}.fasta -I ${ALIGNED_DIR}/deduplicated_${SRA}.bam -O ${VARIANT_DIR}/raw_variants.vcf

# Filter Variants
echo FILTERING VARIANTS...
gatk VariantFiltration -R $RAW_DIR/reference_${REF_ID}.fasta -V $VARIANT_DIR/raw_variants.vcf -O $VARIANT_DIR/filtered_variants.vcf --filter-expression "QD < 2.0 || FS > 60.0" --filter-name FILTER
# filtered variants same as raw -> because our disease is an SNP


# -------------------------------PART 7: FUNCTIONAL ANNOTATION -----------------------------------------------

# Re-download the reference sequence as GenBank format, containing indexed genes, transcripts, exons, CDS, UTRs, etc. 
# We are using our chosen reference seq as the database instead of using snpEff's provided database.
echo Downloading reference GenBank file for snpEff...
efetch -db nucleotide -id $REF_ID -format genbank > $SNPEFF_DATA_DIR/genes.gbk
echo Downloaded GenBank file for snpEff!


# snpEff need the following:
# 1. A reference genome sequence (FASTA) 
# 2. A gene annotation file (GenBank or GTF/GFF3)
# 3. A properly configured snpEff.config file

#   The contig file should include 3 lines:
#       1. reference_db.genome : reference_db	    This gives your custom genome a name (reference_db) to be used later like snpEff -v reference_db your_variants.vcf.
#       2. reference_db.fa : .../reference.fasta	This tells snpEff where your FASTA file is (used to build internal genome structure).
#       3. reference_db.genbank : .../genes.gbk 	This points to a GenBank file that contains gene annotations (exons, CDS, etc.) for the same FASTA file.

# Create snpEff condiguration file; clean up later by making new directory for config files 
echo CREATING CUSTOM SNPEFF CONFIGURATION FILE...
# <<EIF tells cat to start reading the next lines as block of text, until hitting a line that says "EOF". We can actually use any word instead of EOF.
cat <<EOF > $SNPEFF_DIR/snpEff.config       # Write the next few lines to the .config file
# Custom snpEff config for reference_db
reference_db.genome : reference_db
reference_db.fa : $(readlink -f $RAW_DIR/reference_${REF_ID}.fasta)     
reference_db.genbank : $(readlink -f $SNPEFF_DATA_DIR/genes.gbk)
EOF
# readlink is a command that shows where a symbolic link (symlink) points to. 
# The -f option means: follow every symlink in the path and return the absolute path to the final destination file.
# Example: echo $(readlink -f $RAW_DIR/reference_${REF_ID}.fasta)
# Output: /home/user/BINF5002_project/rawdata/reference_NG_059281.fasta


# Build snpEff database
echo Building snpEff database...
snpEff build -c $SNPEFF_DIR/snpEff.config -genbank -v -noCheckProtein reference_db
echo Built snpEff database!

#Exporting snpEff database (this is optional, we're just dumping the contents of the snpEff file into a text file for debugging)
echo Exporting snpEff database...
snpEff dump -c $SNPEFF_DIR/snpEff.config reference_db > $SNPEFF_DIR/snpEff_reference_db.txt 
echo Exported snpEff database!


echo Annotating variants with snpEff...
snpEff -c $SNPEFF_DIR/snpEff.config -stats $SNPEFF_DIR/snpEff.html reference_db $VARIANT_DIR/filtered_variants.vcf > $ANNOTATED_DIR/annotated_variants.vcf
echo Annotated variants with snpEff!




# --------------------------------- SUPPLEMENT ---------------------------------------------
# Use tree to check for correct files and directories.
tree $RESULTS_DIR


# # First do QC check.
# echo Performing QC check on raw FASTQ file...
# fastqc "${RAW_DIR}/${SRA}.fastq" -o $QC_REPORT
# echo QC check completed.

