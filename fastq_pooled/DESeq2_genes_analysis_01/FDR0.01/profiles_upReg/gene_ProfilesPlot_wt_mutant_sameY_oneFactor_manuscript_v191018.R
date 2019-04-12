#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage of REC8_HA_Rep1 and other
# chromatin marks around target and random loci

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")

matDir <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/DESeq2_genes_analysis_01/FDR0.01/profiles_upReg/matrices/"
plotDir <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/DESeq2_genes_analysis_01/FDR0.01/profiles_upReg/plots/"

libNames <- c(
              "REC8_HA_Rep1",
              "kss_REC8_HA_Rep1_ChIP_kss_REC8_HA_Rep1_input",

              "SPO11_1_oligos_RPI1",
              "kss_SPO11_1_oligos_RPI34",

              "CHGmeth",
              "kss_CHGmeth",

              "H3K9me2",
              "kss_H3K9me2",
              
              "SPO11_ChIP4",
              "MNase",

              "PolIV_Rep2",
              "PolV",

              "WT_RNAseq_Chris_Rep1",
              "kss_RNAseq_Chris_Rep1",
              "WT_RNAseq_Chris_Rep2",
              "kss_RNAseq_Chris_Rep2"
             )

# Define column mean coverage outfile (mean profiles)
outDFCM <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir, libNames[[x]],
              "_norm_cov_upReg_genes_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir, libNames[[x]],
              "_norm_cov_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Read in target and ranLoc mean coverage profiles
target_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM[[x]][[1]])
})
ranLoc_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM[[x]][[2]])
})

# Redefine names for use in plots
libNames <- c(
              "Wild type",
              "kyp suvh5 suvh6",

              "Wild type",
              "kyp suvh5 suvh6",

              "Wild type",
              "kyp suvh5 suvh6",

              "Wild type",
              "kyp suvh5 suvh6",

              "Wild type",
              "Wild type",

              "Wild type",
              "Wild type",

              "Wild type",
              "kyp suvh5 suvh6",
              "Wild type",
              "kyp suvh5 suvh6"
             )

YlabelNames <- c(
                 "REC8-HA",
                 "SPO11-1-oligos",
                 "CHG methylation",
                 "H3K9me2",
                 "SPO11-1 ChIP",
                 "MNase",
                 "Pol IV",
                 "Pol V",
                 "RNA-seq",
                 "RNA-seq Rep2"
                )

# Plot mean REC8 vs other coverage profiles around target and random loci
pdf(paste0(plotDir, "kss_upReg_FDR0.01_geneProfiles_wt_mutant_RNAseq_REC8_SPO11_1_oligos_CHGmeth_winSize20_manuscript_v191018.pdf"), height = 10, width = 6)
par(mfrow = c(4, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
xplot <- seq(1, length(target_covDat[[1]][,1]), by = 1)

# RNA-seq rep1
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat[[13]][,1],
    mutantDat1  = target_covDat[[14]][,1],
    ranDat1 = ranLoc_covDat[[13]][,1],
    mutantRanDat1 = ranLoc_covDat[[14]][,1],
    Ylabel1 = YlabelNames[9],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topleft",
    legendLabs = c(libNames[13], libNames[14]),
    wtCol =  "blue", mutantCol = "red")
# REC8
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    mutantDat1  = target_covDat[[2]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    mutantRanDat1 = ranLoc_covDat[[2]][,1],
    Ylabel1 = YlabelNames[1],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topleft",
    legendLabs = c(libNames[1], libNames[2]),
    wtCol =  "blue", mutantCol = "red")
# SPO11-1-oligos
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat[[3]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    ranDat1 = ranLoc_covDat[[3]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    Ylabel1 = YlabelNames[2],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "bottomright",
    legendLabs = c(libNames[3], libNames[4]),
    wtCol =  "blue", mutantCol = "red")
# CHGmeth
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat[[5]][,1],
    mutantDat1  = target_covDat[[6]][,1],
    ranDat1 = ranLoc_covDat[[5]][,1],
    mutantRanDat1 = ranLoc_covDat[[6]][,1],
    Ylabel1 = YlabelNames[3],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topleft",
    legendLabs = c(libNames[5], libNames[6]),
    wtCol =  "blue", mutantCol = "red")
dev.off()

