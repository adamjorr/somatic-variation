CONTIGS=scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 \
scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11

BAM=../../alignment.recal.bam
REF=../../e_mel_3/e_mel_3.fa

ALLCALLS=../allcalls.vcf.gz
DENOVOS=../denovos.vcf.gz

BCFTOOLS=bcftools
BEDTOOLS=bedtools
RSCRIPT=Rscript --vanilla

#####################################################################

default: all

all: deduped.vcf.gz goodsites.vcf figures

.PHONY: default all

.DELETE_ON_ERROR:

.SECONDARY:

#####################################################################

# Construct a genome description file for bedtools
genome.txt: $(DENOVOS)
	$(BCFTOOLS) view -h $< | grep '^##contig' | \
		sed -e 's|##contig=<ID=||' -e 's|,length=|\t|' -e 's|>$$||' \
		> $@

# Remove denovos clustered withing 1000 nt of one another
deduped.bed: $(DENOVOS)
	$(BEDTOOLS) merge -c 1 -o count -d 1000 -i $< | awk '$$4==1' > $@

deduped.vcf.gz: $(DENOVOS) deduped.bed
	$(BCFTOOLS) view -o $@ -O z -T deduped.bed $(DENOVOS)
	$(BCFTOOLS) $@

# Generate regions around deduped denovos
deduped.slop.bed: deduped.bed genome.txt
	bedtools slop -b 500 -g genome.txt -i $< > $@

# Extract reads from near denovos
bam/deduped_regions.bam: $(BAM) deduped.slop.bed
	samtools view -O BAM -o $@ -L deduped.slop.bed $<
	samtools index $@

hplen.tsv.gz: ../phasing/phased.vcf.gz
	${RSCRIPT} hplen.R
	bgzip $(basename $@)
	tabix -S 1 -s 1 -b 2 -e 2 $@

#####################################################################
# Select denovos that are of high-quality

goodsites.vcf: deduped.vcf.gz hplen.tsv.gz
	$(BCFTOOLS) annotate -a hplen.tsv.gz -h hplen_header.txt -c CHROM,POS,HPLEN $< | \
	$(BCFTOOLS) view -o $@ -i 'LLS >= -5 & DNP >= 0.99999 & HPLEN >= 500'

# generate region around good denovos
goodsites.slop.bed: goodsites.vcf genome.txt
	$(BEDTOOLS) merge -i $< | $(BEDTOOLS) slop -b 100 -g genome.txt > $@

goodsites.context.txt: goodsites.vcf genome.txt
	$(BEDTOOLS) merge -i $< | $(BEDTOOLS) slop -b 1 -g genome.txt \
		| $(BEDTOOLS) getfasta -tab -fi $(REF) -bed - > $@

#####################################################################
# FIGURES

figures: figures/fig-em-hist.pdf figures/fig-em-xy.pdf figures/fig-em-tree.pdf

.PHONY: figures

figures/fig-%.pdf: figures/fig-%.R goodsites.vcf
	$(RSCRIPT) --vanilla $<

