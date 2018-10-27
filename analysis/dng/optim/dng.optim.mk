CONTIGS=scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 \
scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11

BAM=../../alignment.recal.bam
REF=../../e_mel_3/e_mel_3.fa

BCFTOOLS=bcftools

BEDS=$(addsuffix .bed,$(CONTIGS))
VCFS=$(addsuffix .bcftools.vcf.gz,$(CONTIGS))

##########################################################

default: all

all: bed vcf

.PHONY: all default bed vcf

.DELETE_ON_ERROR:

.SECONDARY:

#########################################################

bed: $(BEDS)

$(BEDS): %.bed : e_mel_3_degenerate-fixed.bed
	grep -e '^$*\b' $< > $@

vcf: $(VCFS)

%.bcftools.vcf.gz: %.bed
	bcftools mpileup -d 10000 -a 'AD,DP' -r $* -T $< -q 3 -Q 13 -O u -f $(REF) $(BAM) | \
	bcftools call -m -A -p 0 -P 0.025 -O z -o $@

loglike.jobs:
	printf 'dng loglike --ped=sampleM.ped --input=%s\n' $(VCFS) > $@

optimized_params.txt: loglike.jobs
	Rscript --vanilla optimize_parameters.R $< | tee $@
