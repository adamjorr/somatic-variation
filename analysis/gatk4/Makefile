SCRIPTDIR=../../scripts

PLOTTREE=$(SCRIPTDIR)/plot_tree.R
FILTWITHREPLICATES=$(SCRIPTDIR)/filt_with_replicates.pl
VCFTOFA=$(SCRIPTDIR)/vcf2fa.sh

%.vcf.gz: %.vcf
	bgzip -c $< >$@

%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

var-calls.vcf: ../alignment.recal.bam e_mel_3/e_mel_3.fa e_mel_3/e_mel_3.dict
	mkdir -p calls; parallel --halt 2 gatk HaplotypeCaller --heterozygosity 0.025 -R $(word 2, $^) -I $< -ploidy 2 -L scaffold_{} -O calls/{}_$@ ::: $(shell seq 1 11); picard SortVcf $(foreach num, $(shell seq 1 11), I=calls/$(num)_$@) O=$@ SEQUENCE_DICTIONARY=$(word 3, $^)

depth-and-het-filter.vcf: var-calls.vcf
	bcftools filter -g 50 -i 'DP <= 500 && ExcessHet <=40' $< | bcftools view -m2 -M2 -v snps -o $@

repeat-filter.vcf: depth-and-het-filter.vcf ../liftover/e_mel_3_repeatmask.bed
	vcftools --vcf $< --exclude-bed $(word 2, $^) --remove-filtered-all --recode --recode-INFO-all --stdout >$@

#using diploidify -v has the effect of removing nonvariable and triallelic sites
replicate-filter-only-variable.vcf: repeat-filter.vcf
	cat $< | $(FILTWITHREPLICATES) -s -g 3 | diploidify.py -t vcf -v >$@

replicate-filter-only-variable.fa: replicate-filter-only-variable.vcf
	cat $< | $(VCFTOFA) | diploidify.py -v >$@

raxmltmp := $(shell mktemp -d --tmpdir=$(CURDIR) raxml_XXXXXX)
%.nwk: %.fa
	raxml -T 8 -f a -s $< -n nwk -m ASC_GTRGAMMA -w $(raxmltmp) --asc-corr=lewis -p 12345 -x 12345 -# 100 && cat $(raxmltmp)/RAxML_bestTree.nwk >$@

%.pdf: %.nwk
	${PLOTTREE} $< $@
