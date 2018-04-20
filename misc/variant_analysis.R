library(VariantAnnotation)
library(IRanges)
library(GenomicRanges)
library(GenomicFeatures)
library(tidyverse)


# genome <- toGRanges("chromosomes.txt")
# bands <- toGRanges("genes.txt")
# withrepeats <- read_tsv("positions_with_repeats.txt", col_names=c("chr","pos"))
# norepeats <- read_tsv("positions_no_repeats.txt", col_names=c("chr","pos"))
# genes <- read_tsv("genes.txt") %>%
#   group_by(chr)
# 
# in_range_ <- function(chrm, pos, genes){
#   rnges <- genes %>% filter(chr == chrm)
#   btwn <- function(pos, start, end){between(pos,start,end)}
#   result <- map2(rnges$start,rnges$end, ~ btwn(pos,.x,.y))
#   return(TRUE %in% result)
# }
# in_range <- function(chrm,pos,genes){
#   return(map2(chrm,pos,~ in_range_(.x,.y,genes)))
# }
# 
# variant_table <- norepeats %>%
#   mutate(repeatregion = FALSE) %>%
#   right_join(withrepeats) %>%
#   mutate(repeatregion = is.na(.$repeatregion)) %>%
#   mutate(ingene = unlist(in_range(chr,pos,genes)))

vcffile <- "passed.vcf"
fastafile <- "e_mel_3.fa"
gfffile <- "e_mel_3_genes.gff3"
repeatsgfffile <- "e_mel_3_repeats.gff3"
# this may break if tidyverse is loaded??
# vcf <- readVcf(vcffile,fastafile)
gff <- GenomicFeatures::makeTxDbFromGFF(gfffile)
locs <- locateVariants(vcf,gff,AllVariants())
fa <- FaFile(fastafile)
indexFa(fastafile)
coding <- predictCoding(vcf,gff,seqSource=fa)

#overlapsAny(locs, genes(gff)) seems to do what i want


# load tidy things here
positions <- as_tibble(start(ranges(coding)))
names(positions) <- "pos"
chrs <- as_tibble(seqnames(coding))
names(chrs) <- "chr"
coding <- chrs %>%
  bind_cols(positions) %>%
  bind_cols(as_tibble(mcols(coding)))

positions <- as_tibble(start(ranges(locs)))
names(positions) <- "pos"
chrs <- as_tibble(seqnames(locs))
names(chrs) <- "chr"

locs <- chrs %>%
  bind_cols(positions) %>%
  bind_cols(as_tibble(mcols(locs)))

type <- locs %>%
  select(chr, pos, LOCATION) %>%
  distinct()

coding <- coding %>%
  select(chr, pos, CONSEQUENCE, REFAA, VARAA)

variant_table <- variant_table %>%
  left_join(type,by = c("chr","pos")) %>%
  left_join(coding,by = c("chr","pos")) 

# variant_table <- variant_table %>%
#   mutate(seq = unlist(map2(.$chr,.$pos, ~ as.character(getSeq(fa, GRanges(seqnames = .x, IRanges(start = .y - 1000, end = .y + 1000)))))))

gts <- as.tibble(geno(genotypeCodesToNucleotides(vcf))$GT)

refalt <- as_tibble(mcols(rowRanges(vcf))) %>%
  select(REF, ALT) %>%
  mutate(ALT = unlist(map(.$ALT, ~ as.character(.)))) %>%
  bind_cols(as_tibble(as.data.frame(seqnames(rowRanges(vcf))))) %>%
  rename(chr = value) %>%
  bind_cols(as_tibble(start(ranges(rowRanges(vcf))))) %>%
  rename(pos = value) %>%
  select(chr, pos, REF, ALT) %>%
  bind_cols(gts) %>%
  mutate(upstreamkb = unlist(map2(.$chr,.$pos, ~as.character(getSeq(fa,GRanges(seqnames = .x, IRanges(start = .y - 1000, end = .y - 1))))))) %>%
  mutate(downstreamkb = unlist(map2(.$chr,.$pos, ~as.character(getSeq(fa,GRanges(seqnames = .x, IRanges(start = .y + 1 , end = .y + 1000))))))) %>%
  select(chr:ALT, upstreamkb, M1a:M8c, downstreamkb)

variant_table <- variant_table %>%
  left_join(refalt, by=c("chr","pos"))

write_tsv(variant_table,"variants.tsv")
