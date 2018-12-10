INDUCEDMUTRESULTS=../../scripts/spiked_mutation_results.py

.SECONDARY:

#this is from ../2018-05-08_induced_mutations/alignment.bam
#recal.bam: recal.table dedup.bam ../e_mel_3/e_mel_3.fa
#	gatk ApplyBQSR -I $(word 2, $^) -R $(word 3, $^) --bqsr-recal-file $< -O $@

filter/goodsites.vcf: denovos.vcf.gz
	$(MAKE) -C $(dir $@) $(notdir $@)

results.txt: ../alignment.bam ../alignment.bam.bai ../mut_files/mutfile.txt filter/goodsites.vcf inrepeats.txt
	$(INDUCEDMUTRESULTS) -s $< -v $(word 4, $^) -r $(word 5, $^) --dng > $@

mut_files/mutations.tab.gz: ../mut_files/mutations.tab
	bgzip -c <$< >$@

mut_files/mutations.tab.gz.tbi: ../mut_files/mutations.tab.gz
	tabix -b2 -e2 $<

inrepeats.txt: ../mut_files/mutations.tab.gz ../e_mel_3_repeatmask.bed
	tabix -R $(word 2, $^) $< >$@

include dng_caller.mk
