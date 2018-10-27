CONTIGS=scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 \
scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11

BAM=../2018-05-08_induced_mutations/recal.bam
REF=../e_mel_3/e_mel_3.fa

GATK=gatk
DNGCALL=dng call

##########################################################

default: all

all: gatk call \
	allgatk.vcf.gz allgatk.vcf.gz.csi \
	allcalls.vcf.gz allcalls.vcf.gz.csi \
	denovos.vcf.gz denovos.vcf.gz.csi

.PHONY: all default

.DELETE_ON_ERROR:

.SECONDARY:

#########################################################
# Helper Rules

%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

#########################################################
# Run Haplotypecaller to generate a base callset for our data

.PHONY: gatk

gatk: $(addprefix gatk/,$(addsuffix .haplotypecaller.vcf.gz,$(CONTIGS)))

gatk/%.haplotypecaller.vcf.gz: $(BAM)
	$(GATK) HaplotypeCaller -O $@ -R $(REF) -I $(BAM) -L $* \
	--verbosity=ERROR \
	--standard-min-confidence-threshold-for-calling=0 \
	--heterozygosity=0.025 \
	--indel-heterozygosity=0.003125 \
	--base-quality-score-threshold=13 \
	--minimum-mapping-quality=3 \
	--read-filter=PrimaryLineReadFilter

##########################################################
# Run dng-call to call variants and all de novos
#
# Parameters were estimated from bcftools mpileup+call
# output using 3-fold degenerate sites in the ./optim directory

.PHONY: call

call: $(addprefix call/,$(addsuffix .call.vcf.gz,$(CONTIGS)))

call/%.call.vcf.gz call/%.call.vcf.gz.csi: gatk/%.haplotypecaller.vcf.gz
	$(DNGCALL) --all --input=$< --output=$@ \
		--ped=sampleM.ped \
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

allgatk.vcf.gz: $(addprefix gatk/,$(addsuffix .haplotypecaller.vcf.gz,$(CONTIGS)))
	bcftools concat -O z -o $@ --threads 4 $^
	bcftools index $@

temp/allcalls.vcf.gz: $(addprefix call/,$(addsuffix .call.vcf.gz,$(CONTIGS)))
	bcftools concat -O z -o $@ --threads 4 $^
	bcftools index $@

.INTERMEDIATE: temp/allcalls.vcf.gz

allcalls.vcf.gz: allgatk.vcf.gz temp/allcalls.vcf.gz temp/allcalls.vcf.gz.csi
	bcftools annotate -O z -o $@ --threads 4 -a allgatk.vcf.gz -c '+INFO' temp/allcalls.vcf.gz
	bcftools index $@

denovos.vcf.gz: allcalls.vcf.gz
	bcftools view -O z -o $@ -i 'DENOVO=1' $<
	bcftools index $@
