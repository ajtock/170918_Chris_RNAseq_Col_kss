#!/applications/R/R-3.5.0/bin/Rscript

# Use Mann-Whitney-Wilcoxon (U) tests to determine whether significant
# differences exist between genotypes with regard to library-size-normalised
# ChIP-seq coverage within windows along upregulated features and
# non-differentially expressed features

# Usage from within parent directory containing "matrices" directory
# of window coverage files
# ./Utest_coverage_per_win_features.R REC8_HA_Rep2 kss_REC8_HA_Rep1 'wt REC8-HA_Rep2,kss REC8-HA_Rep1' 'red,red4' 2000 2kb 20 20bp 0.1 Copia_LTR

#geno1Name <- "REC8_HA_Rep2"
#geno2Name <- "kss_REC8_HA_Rep1"
#libNamesPlot <- unlist(strsplit("wt REC8-HA_Rep2,kss REC8-HA_Rep1",
#                                split = ","))
#colours <- unlist(strsplit("red,red4",
#                           split = ","))
#flankSize <- 2000
#flankName <- "2kb"
#winSize <- 20
#winName <- "20bp"
#FDR <- 0.1
#FDRname <- paste0("FDR", as.character(FDR))
#featureName <- "Copia_LTR"

args <- commandArgs(trailingOnly = TRUE)
geno1Name <- args[1]
geno2Name <- args[2]
libNamesPlot <- unlist(strsplit(args[3],
                                split = ","))
colours <- unlist(strsplit(args[4],
                           split = ","))
flankSize <- as.numeric(args[5])
flankName <- as.character(args[6])
winSize <- as.numeric(args[7])
winName <- as.character(args[8])
FDR <- as.numeric(args[9])
FDRname <- paste0("FDR", as.character(FDR))
featureName <- args[10]

library(MASS) 

matDir <- "matrices/"
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))

# Load matrices in which each row corresponds to a feature or random locus
# and each column corresponds to a window within that locus
geno1_matList <- list(read.table(paste0(matDir,
                                        geno1Name, "_norm_cov_",
                                        featureName, "_TEs_smoothed_target_and_",
                                        flankName, "_flank_dataframe.txt"),
                                 header = T),
                      read.table(paste0(matDir,
                                        geno1Name, "_norm_cov_",
                                        featureName, "_ranLoc_smoothed_target_and_",
                                        flankName, "_flank_dataframe.txt"),
                                 header = T))
geno2_matList <- list(read.table(paste0(matDir,
                                        geno2Name, "_norm_cov_",
                                        featureName, "_TEs_smoothed_target_and_",
                                        flankName, "_flank_dataframe.txt"),
                                 header = T),
                      read.table(paste0(matDir,
                                        geno2Name, "_norm_cov_",
                                        featureName, "_ranLoc_smoothed_target_and_",
                                        flankName, "_flank_dataframe.txt"),
                                 header = T))

# Define and create new directories and subdirectories to contain results
sizeDir <- paste0(flankName, "_flank_", winName, "_win/")
outDir <- paste0(sizeDir, "Utests/")
plotDir <- paste0(outDir, FDRname, "/")
system(paste0("[ -d ", sizeDir, " ] || mkdir ", sizeDir))
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# For each genotype, calculate mean and variances of coverage
# in each window along features (list element [[1]]) and
# random loci (list element [[2]])
geno1_means <- list(as.vector(apply(X = geno1_matList[[1]],
                                    MARGIN = 2,
                                    FUN = mean)),
                    as.vector(apply(X = geno1_matList[[2]],
                                    MARGIN = 2,
                                    FUN = mean)))
geno2_means <- list(as.vector(apply(X = geno2_matList[[1]],
                                    MARGIN = 2,
                                    FUN = mean)),
                    as.vector(apply(X = geno2_matList[[2]],
                                    MARGIN = 2,
                                    FUN = mean)))

geno1_vars <- list(as.vector(apply(X = geno1_matList[[1]],
                                   MARGIN = 2,
                                   FUN = var)),
                   as.vector(apply(X = geno1_matList[[2]],
                                   MARGIN = 2,
                                   FUN = var)))
geno2_vars <- list(as.vector(apply(X = geno2_matList[[1]],
                                   MARGIN = 2,
                                   FUN = var)),
                   as.vector(apply(X = geno2_matList[[2]],
                                   MARGIN = 2,
                                   FUN = var)))

# Plot means vs variances as a quick check for overdispersion
# data look overdispersed (greater variances than means)
pdf(paste0(outDir, geno1Name, "_", geno2Name,
           "_", featureName, "_upreg_TEs_and_nonDE_TEs_coverage_means_vs_variances_",
           flankName, "_flank_", winName, "_win.pdf"),
    height = 10, width = 10)
par(mfrow = c(2, 2))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = geno1_means[[1]], y = geno1_vars[[1]], pch = 19,
     xlab = "Mean", ylab = "Variance",
     xlim = c(min(geno1_means[[1]], geno1_vars[[1]]),
              max(geno1_means[[1]], geno1_vars[[1]])),
     ylim = c(min(geno1_means[[1]], geno1_vars[[1]]),
              max(geno2_means[[1]], geno1_vars[[1]])),
     main = paste0(geno1Name, " coverage around ", featureName, " upreg TEs"),
     cex.main = 0.8)
plot(x = geno2_means[[1]], y = geno2_vars[[1]], pch = 19,
     xlab = "Mean", ylab = "Variance",
     xlim = c(min(geno2_means[[1]], geno2_vars[[1]]),
              max(geno2_means[[1]], geno2_vars[[1]])),
     ylim = c(min(geno2_means[[1]], geno2_vars[[1]]),
              max(geno2_means[[1]], geno2_vars[[1]])),
     main = paste0(geno2Name, " coverage around ", featureName, " upreg TEs"),
     cex.main = 0.8)
plot(x = geno1_means[[2]], y = geno1_vars[[2]], pch = 19,
     xlab = "Mean", ylab = "Variance",
     xlim = c(min(geno1_means[[2]], geno1_vars[[2]]),
              max(geno1_means[[2]], geno1_vars[[2]])),
     ylim = c(min(geno1_means[[2]], geno1_vars[[2]]),
              max(geno1_means[[2]], geno1_vars[[2]])),
     main = paste0(geno1Name, " coverage around ", featureName, " non-DE TEs"),
     cex.main = 0.8)
plot(x = geno2_means[[2]], y = geno2_vars[[2]], pch = 19,
     xlab = "Mean", ylab = "Variance",
     xlim = c(min(geno2_means[[2]], geno2_vars[[2]]),
              max(geno2_means[[2]], geno2_vars[[2]])),
     ylim = c(min(geno2_means[[2]], geno2_vars[[2]]),
              max(geno2_means[[2]], geno2_vars[[2]])),
     main = paste0(geno2Name, " coverage around ", featureName, " non-DE TEs"),
     cex.main = 0.8)
dev.off()

# Create a factor with two levels (geno1Name and geno2Name)
# to be used as the predictor in regression models
genotype <- factor(c(rep(geno1Name, times = dim(geno1_matList[[1]])[1]),
                     rep(geno2Name, times = dim(geno2_matList[[1]])[1])),
                   levels = c(geno1Name, geno2Name))


# For each window, fit various regression models in which genotype is the
# predictor variable and coverage is the response variable
# Determine which model has best fit

# Normal linear regression, equivalent to a t-test:
genotype_normal <- lapply(seq_along(geno1_matList[[1]]), function(x) {
  glm(formula = c(geno1_matList[[1]][[x]], geno2_matList[[1]][[x]]) ~ genotype,
      family = gaussian(link = "identity"))

})
# Get Bayesian Information Criterion (BIC); the smaller, the better the fit
genotype_normal_BIC <- sapply(seq_along(genotype_normal), function(x) {
  if( !is.infinite(BIC(genotype_normal[[x]])) ) {
    BIC(genotype_normal[[x]])
  } else {
    NA
  }
})
# Get Akaike Information Criterion (AIC); the smaller, the better the fit
genotype_normal_AIC <- sapply(seq_along(genotype_normal), function(x) {
  if( !is.infinite(AIC(genotype_normal[[x]], k = 2)) ) {
    AIC(genotype_normal[[x]], k = 2)
  } else {
    NA
  }
})

# Evaluate model goodness-of-fit using chi-squared test based on the
# residual deviance and degrees of freedom
# A P-value > 0.05 indicates that the model fits the data
genotype_normal_chisqPvals <- sapply(seq_along(genotype_normal), function(x) {
  1 - pchisq(summary(genotype_normal[[x]])$deviance,
             summary(genotype_normal[[x]])$df.residual)
})
print(paste0("P-values from chi-squared goodness-of-fit tests evaluating genotype_normal linear regression models for ",
             length(genotype_normal), " genomic windows:"))
print(genotype_normal_chisqPvals)
print(paste0("Sum of P-values from chi-squared goodness-of-fit tests evaluating genotype_normal linear regression models for ",
             length(genotype_normal), " genomic windows:"))
print(sum(genotype_normal_chisqPvals))
#[1] 0
# Normal regression does not fit the data

# Use the model to predict mean coverage for each genotype
# and corresponding standard errors
# The model predicts the correct mean coverage
genotype_normal_predict <- lapply(seq_along(genotype_normal), function(x) {
  cbind(data.frame(genotype = c(geno1Name, geno2Name)),
        mean = predict(genotype_normal[[x]],
                       newdata = data.frame(genotype = c(geno1Name, geno2Name)),
                       type = "response"),
        SE = predict(genotype_normal[[x]],
                     newdata = data.frame(genotype = c(geno1Name, geno2Name)),
                     type = "response",
                     se.fit = TRUE)$se.fit
       )
})
genotype_normal_correct_predictions <- sapply(seq_along(genotype_normal_predict), function(x) {
  c(print(all.equal(genotype_normal_predict[x][[1]]$mean[1],
                    geno1_means[[1]][x])),
    print(all.equal(genotype_normal_predict[x][[1]]$mean[2],
                    geno2_means[[1]][x])))
})
print(paste0(as.character(sum(genotype_normal_correct_predictions)),
             "/", as.character(length(genotype_normal_correct_predictions))))

# Get R2
genotype_normal_R2 <- sapply(seq_along(genotype_normal), function(x) {
  1 - (genotype_normal[[x]]$deviance/genotype_normal[[x]]$null.deviance)
})
# Get Cohen's F2 as a measure of effect size for linear regression
genotype_normal_F2 <- sapply(seq_along(genotype_normal), function(x) {
  genotype_normal_R2[x] / (1 - genotype_normal_R2[x])
})


## poisson
#genotype_poisson <- lapply(seq_along(geno1_matList[[1]]), function(x) {
#  glm(formula = c(geno1_matList[[1]][[x]], geno2_matList[[1]][[x]]) ~ genotype,
#      family = poisson(link = "log"))
#})
#Error in eval(family$initialize) :
#  negative values not allowed for the 'Poisson' family


## Negative binomial (to account for overdispersion)
#genotype_negbin <- lapply(seq_along(geno1_matList[[1]]), function(x) {
#  glm.nb(formula = c(geno1_matList[[1]][[x]], geno2_matList[[1]][[x]]) ~ genotype,
#         control = glm.control(maxit = 2000))
#})
#Error in eval(family$initialize) :
#  negative values not allowed for the 'Poisson' family

# Alternative would be to use read counts and include library size as an offset; see
# http://www.talkstats.com/threads/non-integers-in-negative-binomial-regression.58082/
# https://stats.stackexchange.com/questions/237963/how-to-formulate-the-offset-of-a-glm

# Or perform Mann-Whitney-Wilcoxon (U) tests comparing geno1Name and geno2Name coverage
# in each window
# features
Utests_feature <- lapply(seq_along(geno1_matList[[1]]), function(x) {
  wilcox.test(x = geno1_matList[[1]][[x]],
              y = geno2_matList[[1]][[x]],
              alternative = "two.sided")
})
# "cannot compute exact p-value with ties"
UtestPvals_feature <- sapply(seq_along(Utests_feature), function(x) {
  Utests_feature[[x]]$p.val
})
# Adjust P-values by applying the Benjamini-Hochberg multiple testing correction
UtestAdjPvals_feature <- p.adjust(p = UtestPvals_feature, method = "BH")

# ranLoc
Utests_ranLoc <- lapply(seq_along(geno1_matList[[2]]), function(x) {
  wilcox.test(x = geno1_matList[[2]][[x]],
              y = geno2_matList[[2]][[x]],
              alternative = "two.sided")
})
# "cannot compute exact p-value with ties"
UtestPvals_ranLoc <- sapply(seq_along(Utests_ranLoc), function(x) {
  Utests_ranLoc[[x]]$p.val
})
# Adjust P-values by applying the Benjamini-Hochberg multiple testing correction
UtestAdjPvals_ranLoc <- p.adjust(p = UtestPvals_ranLoc, method = "BH")


# Calculate difference in coverage between geno2 and geno1 as
# a proportion of coverage in geno1
propChange_feature <- sapply(seq_along(geno1_matList[[1]]), function(x) {
  ( as.vector(colMeans(geno2_matList[[1]])[x]) -
    as.vector(colMeans(geno1_matList[[1]])[x]) ) /
    as.vector(colMeans(geno1_matList[[1]])[x])
})
propChange_ranLoc <- sapply(seq_along(geno1_matList[[2]]), function(x) {
  ( as.vector(colMeans(geno2_matList[[2]])[x]) -
    as.vector(colMeans(geno1_matList[[2]])[x]) ) /
    as.vector(colMeans(geno1_matList[[2]])[x])
})

# Calculate log2(geno2/geno1)
log2FC_feature <- sapply(seq_along(geno1_matList[[1]]), function(x) {
  log2( as.vector(colMeans(geno2_matList[[1]])[x]) /
        as.vector(colMeans(geno1_matList[[1]])[x]) )
})
log2FC_ranLoc <- sapply(seq_along(geno1_matList[[2]]), function(x) {
  log2( as.vector(colMeans(geno2_matList[[2]])[x]) /
        as.vector(colMeans(geno1_matList[[2]])[x]) )
})

# Function for plotting average coverage profiles
plotAvgCov <- function(dat1, dat2,
                       ran1, ran2,
                       col1, col2,
                       mainTitle,
                       flankSize, winSize,
                       flankLabL, flankLabR,
                       startLab, endLab,
                       legendLoc, legendLabs) {
  plot(x = 1:length(dat1),
       y = dat1, col = col1,
       type = "l", lwd = 3, ann = F,
       ylim = c(min(dat1, dat2, ran1, ran2),
                max(dat1, dat2, ran1, ran2)),
       xaxt = "n", yaxt = "n")
  lines(x = 1:length(dat2),
        y = dat2, col = col2,
        type = "l", lwd = 3)
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle)
  axis(side = 2, at = pretty(c(dat1, dat2, ran1, ran2)), lwd = 2, cex.axis = 0.8)
  mtext(side = 2, line = 1.75, cex = 1.0,
        text = bquote("Normalized ChIP"))
  axis(side = 1, lwd = 2,
       at = c(1,
              (flankSize/winSize),
              length(dat1) - (flankSize/winSize),
              length(dat1)),
       labels = c("", "", "", ""))
  mtext(side = 1, line = 0.5, cex = 0.8,
        at = c(1,
               (flankSize/winSize),
               length(dat1) - (flankSize/winSize),
               length(dat1)),
        text = c(flankLabL, startLab, endLab, flankLabR))
  abline(v = c((flankSize/winSize), length(dat1)-(flankSize/winSize)), lty = 3, lwd = 2)
  box(lwd = 2)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2),
         text.col = c(col1, col2),
         text.font = c(1, 1),
         ncol = 1, cex = 0.7, lwd = 2, bty = "n")
}

# Function to plot profile of -log10-transformed adjusted P-values and
# diffrence between kss and wt coverage as a proportion of wt coverage
minusLog10AdjPvalEstimatePlot <- function(AdjPvals1, ests1,
                                          AdjPvals2, ests2,
                                          AdjPvalsCol, estsCol,
                                          mainTitle,
                                          flankSize, winSize,
                                          flankLabL, flankLabR,
                                          startLab, endLab) {
  plot(x = 1:length(AdjPvals1),
       y = -log10(AdjPvals1), col = AdjPvalsCol,
       type = "l", lwd = 3, ann = F,
       ylim = c(0,
                pmax(-log10(0.05), max(c(-log10(AdjPvals1), -log10(AdjPvals2)), na.rm = T))),
       xaxt = "n", yaxt = "n")
  mtext(side = 3, line = 0.5, cex = 1.0, text = mainTitle)
  mtext(side = 2, line = 1.75, cex = 1.0, col = AdjPvalsCol,
        text = bquote("-Log"[10]*"(BH-adjusted "*italic("P")*"-value)"))
  axis(side = 2, cex.axis = 0.8, lwd = 2)
  abline(h = -log10(0.1), lty = 5, lwd = 2, col = AdjPvalsCol)
  par(new = T)
  plot(x = 1:length(ests1),
       y = ests1, col = estsCol,
       type = "l", lwd = 3,
       ylim = c(min(c(ests1, ests2)), max(c(ests1, ests2))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = "",
       cex.main = 1)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.2, adj = c(0.5, -2), xpd = NA, srt = -90, col = estsCol,
       labels = bquote("Log"[2] *
                       "(" * italic("kss") * "/wt)"))
  axis(side = 4, cex.axis = 0.8, lwd = 2)
  axis(side = 1, lwd = 2,
       at = c(1,
              (flankSize/winSize),
              length(AdjPvals1) - (flankSize/winSize),
              length(AdjPvals1)),
       labels = c("", "", "", ""))
  mtext(side = 1, line = 0.5, cex = 0.8,
        at = c(1,
               (flankSize/winSize),
               length(AdjPvals1) - (flankSize/winSize),
               length(AdjPvals1)),
        text = c(flankLabL, startLab, endLab, flankLabR))
  abline(v = c((flankSize/winSize), length(AdjPvals1)-(flankSize/winSize)), lty = 3, lwd = 2)
  box(lwd = 2)
}


# Plot
pdf(paste0(plotDir, geno1Name, "_v_", geno2Name,
           "_ChIP_coverage_around_",
           featureName, "_upreg_and_nonDE_TEs_",
           flankName,  "_flank_", winName, "_win.pdf"),
    height = 8, width = 9)
par(mfrow = c(2, 2))
par(mar = c(2.1, 3.3, 2.1, 3.3))
par(mgp = c(3.0, 0.75, 0))
plotAvgCov(dat1 = colMeans(geno1_matList[[1]]),
           dat2 = colMeans(geno2_matList[[1]]),
           ran1 = colMeans(geno1_matList[[2]]),
           ran2 = colMeans(geno2_matList[[2]]),
           col1 = colours[1],
           col2 = colours[2],
           mainTitle = bquote("Upregulated" ~ .(featureName) ~
                              "(" * italic("n") ~ "=" ~
                              .(prettyNum(dim(geno1_matList[[1]])[1],
                                          big.mark = ",", trim = T)) *
                              ")"),
           flankSize = flankSize,
           winSize = winSize,
           flankLabL = paste0("-", as.character(flankSize/1000), " kb"),
           flankLabR = paste0("+", as.character(flankSize/1000), " kb"),
           startLab = "Start",
           endLab = "End",
           legendLoc = "topleft",
           legendLabs = c(bquote(.(libNamesPlot[1])),
                          as.expression(bquote(italic(.(unlist(strsplit(libNamesPlot[2],
                                                                        split = " "))[1])) ~
                                               .(unlist(strsplit(libNamesPlot[2],
                                                                 split = " "))[2])))))
plotAvgCov(dat1 = colMeans(geno1_matList[[2]]),
           dat2 = colMeans(geno2_matList[[2]]),
           ran1 = colMeans(geno1_matList[[1]]),
           ran2 = colMeans(geno2_matList[[1]]),
           col1 = colours[1],
           col2 = colours[2],
           mainTitle = bquote("Non-DE" ~ .(featureName) ~
                              "(" * italic("n") ~ "=" ~
                              .(prettyNum(dim(geno1_matList[[2]])[1],
                                          big.mark = ",", trim = T)) *
                              ")"),
           flankSize = flankSize,
           winSize = winSize,
           flankLabL = paste0("-", as.character(flankSize/1000), " kb"),
           flankLabR = paste0("+", as.character(flankSize/1000), " kb"),
           startLab = "Start",
           endLab = "End",
           legendLoc = "topleft",
           legendLabs = c(bquote(.(libNamesPlot[1])),
                          as.expression(bquote(italic(.(unlist(strsplit(libNamesPlot[2],
                                                                        split = " "))[1])) ~
                                               .(unlist(strsplit(libNamesPlot[2],
                                                                 split = " "))[2])))))
minusLog10AdjPvalEstimatePlot(AdjPvals1 = UtestAdjPvals_feature,
                              AdjPvals2 = UtestAdjPvals_ranLoc,
                              ests1 = log2FC_feature,
                              ests2 = log2FC_ranLoc,
                              AdjPvalsCol = "purple",
                              estsCol = "forestgreen",
                              mainTitle = bquote("Upregulated" ~ .(featureName) ~
                                                 "(" * italic("n") ~ "=" ~
                                                 .(prettyNum(dim(geno1_matList[[1]])[1],
                                                             big.mark = ",", trim = T)) *
                                                 ")"),
                              flankSize = flankSize,
                              winSize = winSize,
                              flankLabL = paste0("-", as.character(flankSize/1000), " kb"),
                              flankLabR = paste0("+", as.character(flankSize/1000), " kb"),
                              startLab = "Start",
                              endLab = "End")
minusLog10AdjPvalEstimatePlot(AdjPvals1 = UtestAdjPvals_ranLoc,
                              AdjPvals2 = UtestAdjPvals_feature,
                              ests1 = log2FC_ranLoc,
                              ests2 = log2FC_feature,
                              AdjPvalsCol = "purple",
                              estsCol = "forestgreen",
                              mainTitle = bquote("Non-DE" ~ .(featureName) ~
                                                 "(" * italic("n") ~ "=" ~
                                                 .(prettyNum(dim(geno1_matList[[2]])[1],
                                                             big.mark = ",", trim = T)) *
                                                 ")"),
                              flankSize = flankSize,
                              winSize = winSize,
                              flankLabL = paste0("-", as.character(flankSize/1000), " kb"),
                              flankLabR = paste0("+", as.character(flankSize/1000), " kb"),
                              startLab = "Start",
                              endLab = "End")
dev.off()
