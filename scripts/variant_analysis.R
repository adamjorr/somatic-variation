#!/usr/bin/env Rscript
library(VariantAnnotation)
library(IRanges)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(ape)

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
tree_str <- "((((M1:698,M2:638)MD:75,M3:796)MC:773,M4:844)MB:372,(((M7:978,M6:928)MG:111,M8:1029)MF:112,M5:1297)ME:358)MA:0;"

#load tidyverse (locateVariants breaks if you load too early)
library(tidyverse)

#fix tree
tree <- read.tree(text=tree_str)
tree$tip.label <- tree$tip.label %>%
  str_replace("^M","")

#tibbles of useful data
gts <- as.tibble(geno(genotypeCodesToNucleotides(vcf))$GT)

#> plyr::count(str_count(unlist(gts[1,]), pattern = unlist(gts[1,1]))) %>% dplyr::top_n(-1, freq) %>% dplyr::select(matches('x'))
#leastcommon_gt <- function(gts){
#  gts %>%
#    unlist() %>%
#    plyr::count() %>%
#    dplyr::top_n(-1, freq) %>%
#    dplyr::select(matches('x'))
#}

#leastcommon <- gts %>%
#  purrr::transpose() %>%
#  map_dfr(leastcommon_gt)

#mutated_node <- gts %>%
#  map2_dfc(leastcommon, ~str_detect(.x, pattern = .y)) %>%
#  purrr::transpose() %>%
#  map(~names(.[unlist(.)])) %>%
#  map(~str_replace_all(., pattern = "^.*/M", replacement = "")) %>%
#  map(~ .[. %in% tree$tip.label]) %>%
#  map( ~keep.tip(tree, .)$node.label[1]) %>%
#  unlist()

leastcommon_gt <- function(gts){
 gts %>%
   unlist() %>%
   plyr::count() %>%
   dplyr::top_n(-1, freq) %>%
   dplyr::select(matches('x')) %>%
   dplyr::slice(1)
}

mutated_node <- info(vcf)$DNL

if (is.null(mutated_node)) {
	#this set of variants is not from denovogear
	leastcommon <- gts %>%
		purrr::transpose() %>%
		map_dfr(leastcommon_gt)

	mutated_node <- gts %>%
		map2_dfc(leastcommon, ~str_detect(.x, pattern = .y)) %>%
		purrr::transpose() %>%
		map(~names(.[unlist(.)])) %>%
		map(~str_replace_all(., pattern = "^M", replacement = "")) %>%
		map(~str_replace_all(., pattern = "[abc]$", replacement = "")) %>%
		map(~ .[. %in% tree$tip.label]) %>%
		map(~keep.tip(tree, .)) %>%
		map(~ifelse(length(.$tip.label) > 1, .$node.label[1], .$tip.label[1])) %>%
		unlist()
}

variant_table <- tibble(chr = as.character(seqnames(vcf)),
                        pos = start(vcf),
                        repeatregion = overlapsAny(vcf, repeatgff),
                        ingene = overlapsAny(vcf, genes(gff)),
                        coding = overlapsAny(vcf, coding),
                        ref = as.character(rowRanges(vcf)$REF),
                        alt = unstrsplit(CharacterList(rowRanges(vcf)$ALT), sep = ",")
                        )

codingdf <- tibble(chr = as.character(seqnames(coding)),
                   pos = start(coding),
                   consequence = coding$CONSEQUENCE,
                   refaa = as.character(coding$REFAA),
                   varaa = as.character(coding$VARAA))

locsdf <- tibble(chr = as.character(seqnames(locs)),
                 pos = start(locs),
                 locs = as.character(locs$LOCATION))

codingdf <- codingdf %>%
	unique()

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
  bind_cols(mutated_node = mutated_node) %>%
  mutate(downstreamkb = unlist(map2(.$chr,.$pos, ~as.character(getSeq(fa,GRanges(seqnames = .x, IRanges(start = .y + 1 , end = .y + 1000)))))))

write_tsv(variant_table,"variants.tsv")
