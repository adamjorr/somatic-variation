#!/usr/bin/env Rscript
library(VariantAnnotation)
library(IRanges)
library(GenomicRanges)
library(GenomicFeatures)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

#file names
vcffile <- args[1]
fastafile <- args[2]
gfffile <- args[3]
repeatsgfffile <- args[4]

#loaded files
vcf <- readVcf(vcffile,fastafile)
gff <- GenomicFeatures::makeTxDbFromGFF(gfffile)
repeatgff <- import.gff3(repeatsgfffile)
locs <- locateVariants(vcf,gff,AllVariants())
fa <- FaFile(fastafile)
coding <- predictCoding(vcf,gff,seqSource=fa)

#tibbles of useful data
gts <- as.tibble(geno(genotypeCodesToNucleotides(vcf))$GT)

variant_table <- tibble(chr = as.character(seqnames(vcf)),
                        pos = start(vcf),
                        repeatregion = overlapsAny(vcf, repeatgff),
                        ingene = overlapsAny(vcf, genes(gff)),
                        coding = overlapsAny(vcf, coding),
                        ref = as.character(rowRanges(vcf)$REF),
                        alt = as.character(unlist(rowRanges(vcf)$ALT))
                        )

codingdf <- tibble(chr = as.character(seqnames(coding)),
                   pos = start(coding),
                   consequence = coding$CONSEQUENCE,
                   refaa = as.character(coding$REFAA),
                   varaa = as.character(coding$VARAA))

locsdf <- tibble(chr = as.character(seqnames(locs)),
                 pos = start(locs),
                 locs = as.character(locs$LOCATION))

locsdf <- locsdf %>%
  unique() %>%
  nest(locs, .key=locs) %>%
  mutate(locs = map(.$locs,unlist)) %>%
  mutate(locs = unlist(map(.$locs, str_c, collapse = ",")))

#join the relevant data into one table and write
variant_table <- variant_table %>%
  left_join(locsdf, on = c("chr","pos")) %>%
  left_join(codingdf, on=c("chr","pos")) %>%
  mutate(upstreamkb = unlist(map2(.$chr,.$pos,~as.character(getSeq(fa,GRanges(seqnames = .x, IRanges(start = .y - 1000, end = .y - 1))))))) %>%
  bind_cols(gts) %>%
  mutate(downstreamkb = unlist(map2(.$chr,.$pos, ~as.character(getSeq(fa,GRanges(seqnames = .x, IRanges(start = .y + 1 , end = .y + 1000)))))))

write_tsv(variant_table,"variants.tsv")
