#!/bin/bash

# Usage:
# bash check_rrna_unmapped.sh input.bam rRNA_index_prefix output_prefix threads

BAM_FILE=$1
RRNA_INDEX=$2       # e.g., path/to/rrna_index (built with bowtie2-build)
OUT_PREFIX=$3
THREADS=$4

echo "[Step 1] Extracting unmapped reads from BAM..."
samtools view -b -f 4 $BAM_FILE > ${OUT_PREFIX}_unmapped.bam

echo "[Step 2] Converting unmapped BAM to FASTQ..."
samtools fastq ${OUT_PREFIX}_unmapped.bam -1 ${OUT_PREFIX}_unmapped_R1.fq.gz -2 ${OUT_PREFIX}_unmapped_R2.fq.gz -0 /dev/null -s /dev/null -n

echo "[Step 3] Mapping unmapped reads to rRNA reference using bowtie2..."
bowtie2 -x $RRNA_INDEX -1 ${OUT_PREFIX}_unmapped_R1.fq.gz -2 ${OUT_PREFIX}_unmapped_R2.fq.gz \
  -p $THREADS --very-sensitive -S ${OUT_PREFIX}_vs_rRNA.sam 2> ${OUT_PREFIX}_vs_rRNA.log

echo "[Step 4] Parsing mapping results..."
samtools view -F 4 ${OUT_PREFIX}_vs_rRNA.sam | wc -l > ${OUT_PREFIX}_rRNA_mapped_count.txt
echo "Done. See ${OUT_PREFIX}_rRNA_mapped_count.txt for number of reads mapped to rRNA."
