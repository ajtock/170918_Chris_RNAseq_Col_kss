#!/bin/bash
# STAR version 2.5.3a
# samtools version 1.3
# R version 3.3.2

# Example usage via condor submission system on hydrogen node7
# csmit -m 14G -c 24 "bash STAR_PE_RNAseq_mapping_pipeline.sh kss_RNAseq_Chris_Rep1"
i=$1

# align reads to reference genome using STAR
STAR --runThreadN 24 \
     --genomeDir /projects/ajt200/TAIR10/STAR_genome_index/ \
     --readFilesIn ./${i}_R1.fastq.gz ./${i}_R2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix ./${i}_ \
     --outFilterMultimapNmax 10 \
     --outMultimapperOrder Random \
     --outFilterMismatchNmax 2 \
     --outSAMattributes All \
     --twopassMode Basic --twopass1readsN -1 \
     --quantMode TranscriptomeSAM GeneCounts

