CONTIGS=scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 \
scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11

BAM=../../alignment.recal.bam
REF=../../e_mel_3/e_mel_3.fa

ALLCALLS=../allcalls.vcf.gz
DENOVOS=../denovos.vcf.gz

BCFTOOLS=bcftools
BEDTOOLS=bedtools
WHATSHAP=whatshap
RSCRIPT=Rscript --vanilla

#####################################################################

default: all

all: whatshap/input whatshap

.PHONY: default all

.DELETE_ON_ERROR:

.SECONDARY:

#####################################################################
# Whatshap Phasing

# extract vcfs near our denovos
input/%.vcf.gz input/%.vcf.gz.csi: $(ALLCALLS)
	$(BCFTOOLS) view -s GL/M -i 'QUAL >= 20 & DENOVO==0 & GT[0]=="0/1"' -r $* $(ALLCALLS) | \
		$(BCFTOOLS) annotate -O z -o $@ -x FORMAT,INFO
	$(BCFTOOLS) index $@

phased/%.vcf.gz phased/%.vcf.gz.csi: input/%.vcf.gz $(BAM)
	$(WHATSHAP) phase --ignore-read-groups $^ | bgzip > $@
	$(BCFTOOLS) index $@

phased.vcf.gz: $(addprefix phased/,$(addsuffix .vcf.gz,$(CONTIGS)))
	$(BCFTOOLS) concat -O z -o $@ $^

%.vcf.gz.csi: %.vcf.gz
	$(BCFTOOLS) index $<

whatshap: phased.vcf.gz phased.vcf.gz.csi $(addprefix phased/,$(addsuffix .vcf.gz,$(CONTIGS)))

whatshap/input: $(addprefix input/,$(addsuffix .vcf.gz,$(CONTIGS)))

.PHONY: whatshap whashap/input
