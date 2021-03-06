#!/bin/bash

csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep2 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_MYC_Rep1 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_kss_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed kss_REC8_HA_Rep1 0.05" & sleep 20;

csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab.bed MNase 0.05" & sleep 20;

csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI1 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI3 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI8 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/log2suvh456SPO11oligoRPI34NakedDNA_norm_allchrs_coverage_coord_tab.bed kss_SPO11_1_oligos_RPI34 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/log2suvh456SPO11oligoRPI35NakedDNA_norm_allchrs_coverage_coord_tab.bed kss_SPO11_1_oligos_RPI35 0.05" & sleep 20;

csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/WT_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Chris_Rep1 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/WT_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Chris_Rep2 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/kss_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed kss_RNAseq_Chris_Rep1 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/kss_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed kss_RNAseq_Chris_Rep2 0.05" & sleep 20;

csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/log2_WT_H3K9me2_ChIP_WT_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed H3K9me2 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K9me2/kss/coverage/log2ChIPinput/log2_kss_H3K9me2_ChIP_kss_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed kss_H3K9me2 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K9me2/cmt3/coverage/log2ChIPinput/log2_cmt3_H3K9me2_ChIP_cmt3_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed cmt3_H3K9me2 0.05" & sleep 20;

csmit -m 3G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed CGmeth 0.05" & sleep 20;
csmit -m 3G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed CHGmeth 0.05" & sleep 20;
csmit -m 3G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed CHHmeth 0.05" & sleep 20;

csmit -m 3G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CG.bed.gr.tab.bed kss_CGmeth 0.05" & sleep 20;
csmit -m 3G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHG.bed.gr.tab.bed kss_CHGmeth 0.05" & sleep 20;
csmit -m 3G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHH.bed.gr.tab.bed kss_CHHmeth 0.05" & sleep 20;

csmit -m 3G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CG.wig.bed.gr.tab.bed cmt3_CGmeth 0.05" & sleep 20;
csmit -m 3G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHG.wig.bed.gr.tab.bed cmt3_CHGmeth 0.05" & sleep 20;
csmit -m 3G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHH.wig.bed.gr.tab.bed cmt3_CHHmeth 0.05" & sleep 20;

csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_MYC_Rep2_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_MYC_Rep2 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/MSH4/WT/coverage/log2ChIPinput/MSH4_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed MSH4 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11_ChIP/WT/coverage/REC8_MYC_Rep2_input/log2ChIPinput/log2_WT_SPO11_ChIP4_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed SPO11_ChIP4 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11_ChIP/WT/coverage/REC8_MYC_Rep2_input/log2ChIPinput/log2_WT_SPO11_ChIP13_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed SPO11_ChIP13 0.05" & sleep 20;

csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep1_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep1 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep2 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep3_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep3 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me1/coverage/log2ChIPinput/WT_H3K4me1_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me1 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me2/coverage/log2ChIPinput/WT_H3K4me2_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me2 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP12_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP12 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP14_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP14 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP15_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP15 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me1/coverage/log2ChIPinput/WT_H3K27me1_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K27me1 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me3/coverage/log2ChIPinput/WT_H3K27me3_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K27me3 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2A_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2A 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AW_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AW 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AX_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AX 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AZ_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AZ 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/PolIV_Law_Jacobsen_2013_Nature/coverage/log2ChIPinput/PolIV_Rep2_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed PolIV_Rep2 0.05" & sleep 20;
csmit -m 10G -c 1 "./gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/PolV_Liu_Jacobsen_2018_NatPlants/coverage/log2ChIPinput/PolV_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed PolV 0.05" & sleep 20;
wait
