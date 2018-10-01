#!/usr/bin/env Rscript
library(karyoploteR)
library(VariantAnnotation)
library(IRanges)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
vcffile <- args[1]
fastafile <- args[2]

vcf <- readVcf(vcffile,fastafile)
fa <- FaFile(fastafile)

genome <- head(scanFaIndex(fa), n = 11)
chrs <- as.character(seqnames(vcf))
pos <- start(vcf)

pdf("chromosome_plot.pdf",10,8)
kp <- plotKaryotype(genome = genome, cex=.75)
kpPoints(kp, chr=chrs, x=pos, y=.2)
dev.off()
