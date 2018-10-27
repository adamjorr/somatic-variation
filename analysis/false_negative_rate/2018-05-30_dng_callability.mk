GATKCALLER=~/bin/zypy/gatkcaller.sh
MUTBAMS=$(shell find -L mutated_bams -name "*.bam" | sort)

repfiltered_only_mutated.vcf.gz: replicate_filtered.vcf.gz replicate_filtered.vcf.gz.csi mut_files/mutations.tab
	bcftools view -R $(word 3, $^) $< -O z -o $@

replicate_filtered.vcf.gz: replicate_filtered.vcf
	bgzip -c $< > $@

replicate_filtered.vcf.gz.csi: replicate_filtered.vcf.gz
	bcftools index $<

#from ../2018-05-08_induced_mutations/
#mut_files/mutations.tab: mut_files/mutfile.txt
#	cut -d' ' -f1,2 $< | tr ' ' '\t' > $@

replicate_filtered.vcf: filtered_bed_excluded.vcf
	cat $< | filt_with_replicates.pl -s -g 3 > $@

filtered_bed_excluded.vcf: filtered.vcf.gz ../liftover/e_mel_3_repeatmask.bed
	vcftools --gzvcf $< --exclude-bed $(word 2, $^) --remove-filtered-all --recode --recode-INFO-all --stdout >$@

filtered.vcf.gz: var-calls-first11.vcf.gz
	bcftools filter -g 50 -i 'DP <= 500 && ExcessHet <=40' $< | bcftools view -m2 -M2 -v snps -O z -o $@

var-calls.vcf.gz: var-calls.vcf
	bgzip -c $< > $@

var-calls.vcf.gz.csi: var-calls.vcf.gz
	bcftools index $<

var-calls-first11.vcf.gz: var-calls.vcf.gz var-calls.vcf.gz.csi
	bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8,scaffold_9,scaffold_10,scaffold_11 $< -O z -o $@

#var-calls.vcf: ../e_mel_3/e_mel_3.fa fixed_alignment.bam
#	${GATKCALLER} -d ./tmp/ -o $@ -r $< -i $(word 2, $^)

.PRECIOUS: fixed_alignment.bam dedup.bam firstcalls.vcf var-calls.vcf.gz ./intervals

#dedup.bam: fixed_alignment.bam
#	picard MarkDuplicates INPUT=$< OUTPUT=$@ METRICS_FILE=$@.metrics.txt CREATE_INDEX=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

#./intervals: ../e_mel_3/e_mel_3.fa
#	rm -rf $@; mkdir -p $@; gatk SplitIntervals -R $< --scatter-count 32 $(foreach num, $(shell seq 1 11), -L scaffold_$(num)) -O $@

#firstcalls.vcf: dedup.bam ../e_mel_3/e_mel_3.fa ../e_mel_3/e_mel_3.dict ./intervals
#	mkdir -p bqsr/intervals; parallel --halt 2 gatk HaplotypeCaller --heterozygosity 0.025 -R $(word 2, $^) -I $< -L {} -stand-call-conf 50 -ploidy 2 -O bqsr/{}_$@ ::: $(wildcard intervals/*.intervals); picard SortVcf $(foreach interval, $(wildcard intervals/*.intervals), I=bqsr/$(interval)_$@) O=$@ SEQUENCE_DICTIONARY=$(word 3, $^)

#this is from ../2018-05-08_induced_mutations/alignment.bam
#recal.bam: recal.table dedup.bam ../e_mel_3/e_mel_3.fa
#	gatk ApplyBQSR -I $(word 2, $^) -R $(word 3, $^) --bqsr-recal-file $< -O $@

#var-calls.vcf: recal.bam ../e_mel_3/e_mel_3.fa ../e_mel_3/e_mel_3.dict
#	mkdir -p calls/intervals; parallel --halt 2 gatk HaplotypeCaller --heterozygosity 0.025 -R $(word 2, $^) -I $< -L {} -ploidy 2 -O calls/{}_$@ ::: $(wildcard intervals/*.intervals); picard SortVcf $(foreach interval, $(wildcard intervals/*.intervals), I=calls/$(interval)_$@) O=$@ SEQUENCE_DICTIONARY=$(word 3, $^)

results.txt: ../2018-05-08_induced_mutations/fixed_alignment.bam ../2018-05-08_induced_mutations/fixed_alignment.bam.bai mut_files/mutfile.txt filter/goodsites.vcf inrepeats.txt
	~/bin/zypy/misc/spiked_mutation_results.py -s $< -v $(word 4, $^) -r $(word 5, $^) --dng > $@

mut_files/mutations.tab.gz: mut_files/mutations.tab
	bgzip -c <$< >$@

mut_files/mutations.tab.gz.tbi: mut_files/mutations.tab.gz
	tabix -b2 -e2 $<

inrepeats.txt: mut_files/mutations.tab.gz ../liftover/e_mel_3_repeatmask.bed
	tabix -R $(word 2, $^) $< >$@

include dng_caller.mk
