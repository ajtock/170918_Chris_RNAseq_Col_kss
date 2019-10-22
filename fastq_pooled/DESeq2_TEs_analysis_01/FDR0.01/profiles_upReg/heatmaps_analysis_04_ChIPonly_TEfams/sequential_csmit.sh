#!/bin/bash

## REC8
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/REC8_MYC_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed REC8_MYC_Rep1 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/REC8_HA_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/coverage/REC8_HA_Rep2_ChIP_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep2 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/kss_REC8_HA_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed kss_REC8_HA_Rep1 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/180830_Chris_lambing_ChIP_rec8HA_kss_rep2/coverage/kss_REC8_HA_Rep2_ChIP_norm_allchrs_coverage_coord_tab.bed kss_REC8_HA_Rep2 0.01" & sleep 10;
#
## SPO11-1-oligos
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/WT_SPO11oligo_RPI1_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI1 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/WT_SPO11oligo_RPI3_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI3 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/WT_SPO11oligo_RPI8_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI8 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/suvh456_SPO11oligo_RPI34_norm_allchrs_coverage_coord_tab.bed kss_SPO11_1_oligos_RPI34 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/suvh456_SPO11oligo_RPI35_norm_allchrs_coverage_coord_tab.bed kss_SPO11_1_oligos_RPI35 0.01" & sleep 10;
#
## H3K9me2
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/WT/coverage/WT_H3K9me2_ChIP_norm_allchrs_coverage_coord_tab.bed H3K9me2 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/kss/coverage/kss_H3K9me2_ChIP_norm_allchrs_coverage_coord_tab.bed kss_H3K9me2 0.01" & sleep 10;
#
# DNA methylation
csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript ./TE_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed mCG 0.01" & sleep 10;
csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript ./TE_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed mCHG 0.01" & sleep 10;
csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript ./TE_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed mCHH 0.01" & sleep 10;

csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript ./TE_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CG.wig.bed.gr.tab.bed kss_mCG 0.01" & sleep 10;
csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript ./TE_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHG.wig.bed.gr.tab.bed kss_mCHG 0.01" & sleep 10;
csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript ./TE_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHH.wig.bed.gr.tab.bed kss_mCHH 0.01" & sleep 10;

csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript ./TE_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CG.wig.bed.gr.tab.bed cmt3_mCG 0.01" & sleep 10;
csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript ./TE_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHG.wig.bed.gr.tab.bed cmt3_mCHG 0.01" & sleep 10;
csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript ./TE_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHH.wig.bed.gr.tab.bed cmt3_mCHH 0.01" & sleep 10;
#
## RNA-seq
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/WT_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Chris_Rep1 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/WT_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Chris_Rep2 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/kss_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed kss_RNAseq_Chris_Rep1 0.01" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/kss_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed kss_RNAseq_Chris_Rep2 0.01" & sleep 10;

wait
