#!/bin/bash

# Define paths
INPUT_DIR="/path/to/input"
OUTPUT_DIR="/path/to/output"
REFERENCE_GENOME="/path/to/bowtie2/reference/genome"
ANNOTATION_FILE="/path/to/annotation.gtf"
QC_DIR="${OUTPUT_DIR}/qc"
TRIMMED_DIR="${OUTPUT_DIR}/trimmed"
ALIGNED_DIR="${OUTPUT_DIR}/aligned"
COUNTS_DIR="${OUTPUT_DIR}/counts"

# Create output directories
mkdir -p $QC_DIR $TRIMMED_DIR $ALIGNED_DIR $COUNTS_DIR

# Step 1: Quality Control (FASTQC)
for FILE in ${INPUT_DIR}/*.fasta; do
    fastqc $FILE -o $QC_DIR
done

# Step 2: Adapter Trimming (Cutadapt)
for FILE in ${INPUT_DIR}/*.fasta; do
    BASENAME=$(basename $FILE .fasta)
    cutadapt -a AGATCGGAAGAGC -o ${TRIMMED_DIR}/${BASENAME}_trimmed.fasta $FILE
done

# Step 3: Alignment (Bowtie2)
for FILE in ${TRIMMED_DIR}/*.fasta; do
    BASENAME=$(basename $FILE _trimmed.fasta)
    bowtie2 -x $REFERENCE_GENOME -U $FILE -S ${ALIGNED_DIR}/${BASENAME}.sam
done

# Step 4: Convert SAM to BAM and Sort
for FILE in ${ALIGNED_DIR}/*.sam; do
    BASENAME=$(basename $FILE .sam)
    samtools view -Sb $FILE | samtools sort -o ${ALIGNED_DIR}/${BASENAME}_sorted.bam
done

# Step 5: Counting with FeatureCounts
featureCounts -t piRNA -T 25 -g piRNA_code  -M -s 1 -F Gtf --minOverlap 15 --fracOverlap 0.95 --fracOverlapFeature 0.95 -a $ANNOTATION_FILE -o ${COUNTS_DIR}/counts.txt ${ALIGNED_DIR}/*_sorted.bam

awk 'BEGIN {OFS="\t"} /^#/ {next} {print $1, $7, $8}' ${COUNTS_DIR}/counts.txt > ${COUNTS_DIR}/count_matrix.txt

