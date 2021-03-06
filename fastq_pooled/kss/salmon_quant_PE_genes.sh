#!/bin/bash
# Quantify gene expression levels using salmon 0.9.1

# Example usage via condor submission system on hydrogen node7
# csmit -m 1G -c 24 "./salmon_quant_PE_genes.sh kss_RNAseq_Chris"

prefix=$1

for i in $prefix"_Rep"{1..2};
do
  samp=`basename ${i}`
  echo "Processing sample ${samp}"
  salmon quant --index /projects/ajt200/TAIR10/salmon_transcriptome_index \
               --libType ISF \
               --mates1 ${i}_R1.fastq.gz \
               --mates2 ${i}_R2.fastq.gz \
               --threads 24 \
               --output ./salmon_quants_genes/${samp}_quant
done
