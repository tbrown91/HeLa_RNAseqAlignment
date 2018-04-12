source("https://bioconductor.org/biocLite.R")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("DESeq2")
library("edgeR")
register(SerialParam())

# Create matrices for gene counts -----------------------------------------

gtffile <- file.path("gencode.v19.chr_patch_hapl_scaff.annotation.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
ebg <- exonsBy(txdb, by="gene")

#Rep1 plus strand
filenames <- list.files(file.path("rep1_plus"))
setwd("rep1_plus")
bamfiles <- BamFileList(filenames)
se1 <- summarizeOverlaps(features=ebg, reads=bamfiles,mode="IntersectionStrict",singleEnd=FALSE,ignore.strand=FALSE,fragments=TRUE )
plus_rep1 = rowSums(assay(se1))
setwd("../")
save(plus_rep1,file="plus_rep1Counts.RData")

#Rep1 minus strand
filenames <- list.files(file.path("rep1_minus"))
setwd("rep1_minus")
bamfiles <- BamFileList(filenames)
se2 <- summarizeOverlaps(features=ebg, reads=bamfiles,mode="IntersectionStrict",singleEnd=FALSE,ignore.strand=FALSE,fragments=TRUE )
minus_rep1 = rowSums(assay(se2))
setwd("../")
save(minus_rep1,file="minus_rep1Counts.RData")

#Rep2 plus strand
filenames <- list.files(file.path("rep2_plus"))
setwd("rep2_plus")
bamfiles <- BamFileList(filenames)
se3 <- summarizeOverlaps(features=ebg, reads=bamfiles,mode="IntersectionStrict",singleEnd=FALSE,ignore.strand=FALSE,fragments=TRUE )
plus_rep2 = rowSums(assay(se3))
setwd("../")
save(plus_rep2,file="plus_rep2Counts.RData")

#Rep2 minus strand
filenames <- list.files(file.path("rep2_minus"))
setwd("rep2_minus")
bamfiles <- BamFileList(filenames)
se4 <- summarizeOverlaps(features=ebg, reads=bamfiles,mode="IntersectionStrict",singleEnd=FALSE,ignore.strand=FALSE,fragments=TRUE )
minus_rep2 = rowSums(assay(se4))
setwd("../")
save(minus_rep2,file="minus_rep2Counts.RData")

rep1 <- plus_rep1 + minus_rep1
rep2 <- plus_rep2 + minus_rep2

r <- cor(rep1,rep2)

HeLa_RNAseqCounts = rep1+rep2
write.table(HeLa_RNAseqCounts,file="HeLa_RNAseqCounts.txt", sep = "\t")

# Plot counts  ------------------------------------------------------------
library("ggplot2")
qplot(log10(rep1),log10(rep2),xlab="Rep1 reads per gene, log10",ylab="Rep2 reads per gene log10")
