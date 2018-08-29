#!/applications/R/R-3.4.0/bin/Rscript
### Convert gene-level counts into format for use with DESeq2

# R version 3.4.0

library(tximport)
print(packageVersion("tximport"))
#[1] ‘1.4.0’

inDir <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled"
outDir <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/DESeq2_genes_analysis_01/"

# Read in table of sample IDs that will be used to specify paths to count files
samples <- read.table(file.path(inDir, "/samples_DESeq2_genes.txt"), header = T)
print(samples)
#  genotype     directory                      sample
#1       WT salmon_quants_genes  WT_RNAseq_Chris_Rep1_quant
#2       WT salmon_quants_genes  WT_RNAseq_Chris_Rep2_quant
#3      kss salmon_quants_genes kss_RNAseq_Chris_Rep1_quant
#4      kss salmon_quants_genes kss_RNAseq_Chris_Rep2_quant

# Specify paths to count files
files <- file.path(inDir, samples$genotype, samples$directory, samples$sample, "quant.sf")
# Set 1:#_samples
names(files) <- paste0("sample", 1:4)
all(file.exists(files))

# Create a dataframe of transcript IDs and corresponding gene IDs
transID <- read.table(files[1], colClasses = c(NA, rep("NULL", 4)), header = T)
tx2gene <- data.frame(cbind(as.vector(transID[,1]), substr(transID[,1], 1, 9)))
colnames(tx2gene) <- c("TXNAME", "GENEID")

# Import transcript-level counts, summarised at gene level
# (reads that map to transcript IDs with a common parent gene ID are pooled)
library(readr)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
print(names(txi))
#[1] "abundance"           "counts"              "length"             
#[4] "countsFromAbundance"

# Import transcript-level counts, summarised at transcript level with "txOut = TRUE"
txi.tx <- tximport(files, type = "salmon", txOut = TRUE, tx2gene = tx2gene)
# Then summarise to gene level
txi.sum <- summarizeToGene(txi.tx, tx2gene)
# These two approaches should produce identical results
print(all.equal(txi$counts, txi.sum$counts))
#[1] TRUE

print(head(txi$counts))
#          sample1 sample2  sample3 sample4
#AT1G01010     305     263   304.00     325
#AT1G01020     293     257   219.00     270
#AT1G01030     447     465   664.00     503
#AT1G01040    2538    1761  1843.00    1467
#AT1G01050     733     600   828.00     481
#AT1G01060   11729    3679 17554.66   16016

save(txi,
     file = paste0(outDir, "tximport.RData"))

