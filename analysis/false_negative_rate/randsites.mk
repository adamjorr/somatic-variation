CONTIGS=scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 \
scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11

BAM=../alignment.recal.bam
REF=../e_mel_3/e_mel_3.fa
BED=sites_to_mutate.bed

BEDS=$(addsuffix .bed,$(CONTIGS))
VCFS=$(addsuffix .bcftools.vcf.gz,$(CONTIGS))

GATK=gatk
DNGCALL=/storage/2017-03-15-eucalyptus/dng/denovogear-3ae70ba8272ca4d6aa96bdb41170ce0739a5aa4a/build/src/dng-call

##########################################################

default: all

all: bed vcf call sites_to_mutate.vcf.gz

.PHONY: all default bed vcf

.DELETE_ON_ERROR:

.SECONDARY:

#########################################################
# Helper Rules

%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

########################################################
# Split BEDS

bed: $(BEDS)

$(BEDS): %.bed : $(BED)
	grep -e '^$*\b' $< > $@

#########################################################
# Make Pileups

vcf: $(VCFS) $(addsuffix .csi,$(VCFS))

%.bcftools.vcf.gz: %.bed
	bcftools mpileup -d 10000 -a 'AD,DP' -r $* -T $< -q 3 -Q 13 -O u -f $(REF) $(BAM) | \
	bcftools call -m -A -p 0 -P 0.025 -O z -o $@

##########################################################
# Run dng-call to call variants and all de novos
#
# Parameters were estimated from bcftools mpileup+call
# output using 3-fold degenerate sites in the ./optim directory

.PHONY: call

CALLS=$(addsuffix .call.vcf.gz,$(CONTIGS))

call: $(CALLS) $(addsuffix .csi,$(CALLS))

%.call.vcf.gz: %.bcftools.vcf.gz
	$(DNGCALL) --all --input=$< --output=$@ \
		--ped=sampleM_star.ped \
		--min-qual=0 \
		--mu=0 \
		--mu-library=0 \
		--mu-somatic=0 \
		--theta=0.0254523129 \
		--lib-error=0.0007602493 \
		--ref-bias-hom=0.9244618851 \
		--ref-bias-het=0.7897195864 \
		--lib-overdisp-hom=0.0074086620 \
		--lib-overdisp-het=0.1914599638 \
		--lib-bias=1

# combine and index all the calls together

sites_to_mutate.vcf.gz sites_to_mutate.vcf.gz.csi: $(CALLS)
	bcftools concat -O z -o $@ --threads 4 $^
	bcftools index $@
