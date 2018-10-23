#!/usr/bin/env Rscript
#http://bioconductor.org/packages/release/bioc/vignettes/signeR/inst/doc/signeR-vignette.html
library(signeR)
library(VariantAnnotation)
library(Rsamtools)

vcffile <- '/storage/2017-03-15-eucalyptus/dng/filter/goodsites.vcf'
fastafile <- '/storage/2017-03-15-eucalyptus/e_mel_3/e_mel_3.fa'

vcf <- readVcf(vcffile, fastafile)
fa <- FaFile(fastafile)

mut <- genCountMatrixFromVcf(fa, vcf)
opp <- genOpportunityFromGenome(fa, scanFaIndex(fa)[1:11])
opp2 <- matrix(opp, nrow = dim(mut)[1], ncol = dim(mut)[2], byrow = TRUE)

signatures <- signeR(M = mut, Opport = opp2, nsig = 1)

#pdf("context.pdf", width = 11, height = 7)
SignPlot(signatures$SignExposures, plot_to_file = TRUE, plots_per_page = 1)
#dev.off()

