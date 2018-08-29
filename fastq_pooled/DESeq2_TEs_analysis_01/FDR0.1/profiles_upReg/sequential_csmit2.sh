#!/bin/bash

csmit -m 10G -c 1 "./TE_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI1 0.1" & sleep 20;

csmit -m 10G -c 1 "./TE_Profiles_commandArgs.R /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/WT_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Chris_Rep2 0.1" & sleep 20;

csmit -m 10G -c 1 "./TE_Profiles_commandArgs.R /projects/ajt200/BAM_masters/PolIV_Law_Jacobsen_2013_Nature/coverage/log2ChIPinput/PolIV_Rep2_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed PolIV_Rep2 0.1" & sleep 20;
wait
