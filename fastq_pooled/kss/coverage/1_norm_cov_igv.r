############################################################################################
# Calculate per base coverage normalised by total genome-wide coverage for each chromosome #
############################################################################################

source("/projects/ajt200/Rscripts/xiaohui/testrun.r")
library(IRanges)
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(doParallel)
# Change number of cores to reflect number of samples you want to process simultaneously
registerDoParallel(cores = 2)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

# Amend as appropriate for your samples
dat.dir <- c("/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/")
# Amend as appropriate for your samples
lib.names1 <- c("kss_RNAseq_Chris_Rep1", "kss_RNAseq_Chris_Rep2")
# Amend as appropriate for your samples
lib.names2 <- c("kss_RNAseq_Chris_Rep1", "kss_RNAseq_Chris_Rep2")
lib.files1 <- paste0(dat.dir, lib.names1, "_STAR_mapped_both_sort.bam")
chr.index <- c(1:5)
# Amend as appropriate for your samples
out.dir <- c("/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/")
cov.files1 <- list(paste0(out.dir, lib.names2[1], "_chr", chr.index, "_unique_both_coverage.txt"),
                   paste0(out.dir, lib.names2[2], "_chr", chr.index, "_unique_both_coverage.txt"))

# Change "foreach(i = 1:2, ... )" to reflect number of samples you want to process simultaneously
foreach(i = 1:2, .combine = 'c') %dopar% {
  print(i)
  index <- i
  test <- chrdat_pair_function(i, lib.files1, cov.files1)
  gc()
}

######################
# Normalise the data #
######################
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chr.lens <- c(30427671,19698289,23459830,18585056,26975502)

# Change "foreach(j = 1:2, ... )" to reflect number of samples you want to process simultaneously
foreach(j = 1:2, .combine = 'c') %dopar% {
 print(lib.names2[j])
 out.files <- paste0(out.dir, lib.names2[j], "_norm_Chr", c(1:5), "_coverage.bedgraph")
 wt.file <- cov.files1[[j]]
 myfiles1 <- read.table(wt.file[1], header=F)
 myfiles2 <- read.table(wt.file[2], header=F)
 myfiles3 <- read.table(wt.file[3], header=F)
 myfiles4 <- read.table(wt.file[4], header=F)
 myfiles5 <- read.table(wt.file[5], header=F)
 totalSum <- colSums(myfiles1)+colSums(myfiles2)+colSums(myfiles3)+colSums(myfiles4)+colSums(myfiles5)
 ##
 myfiles1.new1 <- myfiles1*1e9/totalSum
 myfiles1.new <- cbind(chrs[1], c(0:(chr.lens[1]-1)), c(1:chr.lens[1]), myfiles1.new1)
 write.table(myfiles1.new, file=out.files[1], row.names=F, col.names=F, quote=F)
 ##
 myfiles2.new1 <- myfiles2*1e9/totalSum
 myfiles2.new <- cbind(chrs[2], c(0:(chr.lens[2]-1)), c(1:chr.lens[2]), myfiles2.new1)
 write.table(myfiles2.new, file=out.files[2], row.names=F, col.names=F, quote=F)
 ##
 myfiles3.new1 <- myfiles3*1e9/totalSum
 myfiles3.new <- cbind(chrs[3], c(0:(chr.lens[3]-1)), c(1:chr.lens[3]), myfiles3.new1)
 write.table(myfiles3.new, file=out.files[3], row.names=F, col.names=F, quote=F)
 ##
 myfiles4.new1 <- myfiles4*1e9/totalSum
 myfiles4.new <- cbind(chrs[4], c(0:(chr.lens[4]-1)), c(1:chr.lens[4]), myfiles4.new1)
 write.table(myfiles4.new, file=out.files[4], row.names=F, col.names=F, quote=F)
 ##
 myfiles5.new1 <- myfiles5*1e9/totalSum
 myfiles5.new <- cbind(chrs[5], c(0:(chr.lens[5]-1)), c(1:chr.lens[5]), myfiles5.new1)
 write.table(myfiles5.new, file=out.files[5], row.names=F, col.names=F, quote=F)
}

