CONTIGS=scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 \
scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11

JOBNUMS=$(shell seq 1 100)

BAM=../alignment.recal.bam
REF=../e_mel_3/e_mel_3.fa

DNGDIR=../dng/

DNGCALL=dng call
SHUF_LABELS=../../scripts/label_permutation.R
SHELL=bash
##########################################################

default: all

all: call \
	allcalls.vcf.gz allcalls.vcf.gz.csi \
	denovos.vcf.gz denovos.vcf.gz.csi

.PHONY: all default

.DELETE_ON_ERROR:

.SECONDARY:

#########################################################
# Helper Rules

%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

##########################################################
# Run dng-call to call variants and all de novos
#
# Parameters were estimated from bcftools mpileup+call
# output using 3-fold degenerate sites in the ./optim directory

.PHONY: call

#call: $(addprefix call/,$(addsuffix .call.vcf.gz,$(CONTIGS)))

call: $(addprefix denovos/,$(addsuffix .denovos.vcf.gz, $(JOBNUMS)))

GOODSITES=$(addprefix goodsites/,$(addsuffix .goodsites.vcf.gz, $(JOBNUMS)))

denovos/%.denovos.vcf.gz denovos/%.denovos.vcf.gz.csi: trees/%.ped allgatk.vcf.gz | denovos/
	$(DNGCALL) --input=allgatk.vcf.gz --output=$@ \
		--ped=$< \
		--min-qual=3 \
		--mu=0 \
		--mu-library=0 \
		--mu-somatic=1e-7 \
		--theta=0.0254523129 \
		--lib-error=0.0007602493 \
		--ref-bias-hom=0.9244618851 \
		--ref-bias-het=0.7897195864 \
		--lib-overdisp-hom=0.0074086620 \
		--lib-overdisp-het=0.1914599638 \
		--lib-bias=1
	bcftools index $@

# combine and index all the calls together

#allcalls.vcf.gz: $(addprefix call/,$(addsuffix .call.vcf.gz,$(CONTIGS)))
#	bcftools concat -O z -o $@ --threads 4 $^

allgatk.vcf.gz: $(addprefix $(DNGDIR)/gatk/, $(addsuffix .haplotypecaller.vcf.gz, $(CONTIGS)))
	bcftools concat -O z -o $@ --threads 4 $^

denovos/:
	mkdir -p $@

#denovos/%.denovos.vcf.gz: call/%.call.vcf.gz denovos/
#	bcftools view -O z -o $@ -i 'DENOVO=1' $<

trees/:
	mkdir -p $@

trees/%.ped: | trees/
	echo "##PEDNG v1.0" > $@; echo -e "M\t.\t.\t0" | paste - <( $(SHUF_LABELS) ) >> $@

filter/:
	mkdir -p $@

filter/%.deduped.bed: denovos/%.denovos.vcf.gz | filter/
	bedtools merge -c 1 -o count -d 1000 -i $< | awk '$$4==1' > $@

filter/%.deduped.vcf.gz: denovos/%.denovos.vcf.gz filter/%.deduped.bed
	bcftools view -o $@ -O z -T $(word 2, $^) $<

goodsites/:
	mkdir -p $@

goodsites/%.goodsites.vcf.gz: filter/%.deduped.vcf.gz filter/%.deduped.vcf.gz.csi $(DNGDIR)/filter/hplen.tsv.gz | goodsites/
	bcftools annotate -a $(word 3, $^) -h $(DNGDIR)/filter/hplen_header.txt -c CHROM,POS,HPLEN $< | \
	bcftools view -o $@ -O z -i 'LLS >= -5 & DNP >= 0.99999 & HPLEN >= 500'

num_false_pos.txt: $(GOODSITES)
	rm -f $@; for f in $(GOODSITES) ; do bcftools view -H $$f | wc -l >> $@; done

allgoodsites.bed: $(GOODSITES)
	rm -f $@; for f in $(GOODSITES) ; do bcftools query -f '%CHROM\t%POS0\t%POS\n' $$f >> $@; done

#goodsites_files.txt: $(GOODSITES)
#	for f in goodsites/* ; do echo $$f >> $@; done

num_false_pos_no_originalcalls.txt: $(DNGDIR)/filter/goodsites.vcf $(GOODSITES)
	rm -f $@; for f in $(GOODSITES) ; do bcftools view -H -T ^$< $$f | wc -l >> $@; done

originalcalls_overlaps.vcf.gz: allgoodsites.bed $(DNGDIR)/filter/goodsites.vcf
	bcftools view -T $< -o $@ -O z $(word 2, $^)
