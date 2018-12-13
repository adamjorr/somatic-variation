SCRIPTDIR=../../scripts
INDUCEDMUTRESULTS=$(SCRIPTDIR)/spiked_mutation_results.py

.SECONDARY:

#this is from ../2018-05-08_induced_mutations/alignment.bam
#recal.bam: recal.table dedup.bam ../e_mel_3/e_mel_3.fa
#	gatk ApplyBQSR -I $(word 2, $^) -R $(word 3, $^) --bqsr-recal-file $< -O $@

denovos.vcf.gz:
	$(MAKE) -f ../../dng/Makefile $($@)

filter/goodsites.vcf: denovos.vcf.gz
	$(MAKE) -C $(dir $@) $(notdir $@)

inrepeats.txt: ../mut_files/mutations.tab.gz ../liftover/e_mel_3_repeatmask.bed
	tabix -R $(word 2, $^) $< >$@

results.txt: ../alignment.bam ../alignment.bam.bai ../mut_files/mutfile.txt filter/goodsites.vcf inrepeats.txt
	$(INDUCEDMUTRESULTS) -s $< -v $(word 4, $^) -r $(word 5, $^) --dng > $@

# VCF -> 1k * 14 random positions (fmt CHROM POS VCFREF VCFALT REF ALT) -> split_mutations_by_tree.sh -> induce_mutations.sh
