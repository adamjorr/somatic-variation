GATKCALLER=gatkcaller.sh
PLOTTREE=plot_tree.R

%.vcf.gz: %.vcf
	bgzip -c $< >$@

%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

var-calls.vcf: ../e_mel_3/e_mel_3.fa ../alignment.dedup.bam
	mkdir -p ./tmp/; ${GATKCALLER} -d ./tmp/ -o $@ -r $< -i $(word 2, $^)

first-11-scaffolds.vcf: var-calls.vcf.gz var-calls.vcf.gz.csi
	bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8,scaffold_9,scaffold_10,scaffold_11 -o $@ $<

depth-and-het-filter.vcf: first-11-scaffolds.vcf
	bcftools filter -g 50 -i 'DP <= 500 && ExcessHet <=40' $< | bcftools view -m2 -M2 -v snps -o $@

repeat-filter.vcf: depth-and-het-filter.vcf ../liftover/e_mel_3_repeatmask.bed
	vcftools --vcf $< --exclude-bed $(word 2, $^) --remove-filtered-all --recode --recode-INFO-all --stdout >$@

#using diploidify -v has the effect of removing nonvariable and triallelic sites
replicate-filter-only-variable.vcf: repeat-filter.vcf
	cat $< | filt_with_replicates.pl -s -g 3 | diploidify.py -t vcf -v >$@

replicate-filter-only-variable.fa: replicate-filter-only-variable.vcf
	cat $< | vcf2fa.sh | diploidify.py -v >$@

raxmltmp := $(shell mktemp -d --tmpdir=$(CURDIR) raxml_XXXXXX)
%.nwk: %.fa
	raxml -T 8 -f a -s $< -n nwk -m ASC_GTRGAMMA -w $(raxmltmp) --asc-corr=lewis -p 12345 -x 12345 -# 100 && cat $(raxmltmp)/RAxML_bestTree.nwk >$@

%.pdf: %.nwk
	${PLOTTREE} $< $@
