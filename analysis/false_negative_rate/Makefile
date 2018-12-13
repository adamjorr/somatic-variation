SCRIPTDIR=../scripts
NGMALIGNER=$(SCRIPTDIR)/ngm_aligner.sh
SETUPMUTATIONS=$(SCRIPTDIR)/setup_mutations.sh
SPLITBYTREE=$(SCRIPTDIR)/split_mutations_by_tree.sh
INDUCEMUTATIONS=$(SCRIPTDIR)/induce_mutations.sh

BAM=../alignment.recal.bam
REF=../e_mel_3/e_mel_3.fa

SAMPLENUMBERS=$(shell seq 1 8)
MUTBASE=$(addsuffix .txt, $(addprefix mut_M, $(SAMPLENUMBERS)))
MUTFILES=$(addprefix mut_files/, $(MUTBASE))
SEDSCRIPTS=$(foreach num, 1 2, $(addsuffix .$(num).sed $(addprefix sed_scripts/, $(MUTBASE))))

RAWREADSDIR=../../data/raw
MUTATEDREADSDIR=../../data/mutated


#induce mutations to check false negative rate
# 1k * 14 random positions -> VCF -> sites_to_mutate.txt (fmt CHROM POS VCFREF VCFALT REF ALT) -> split_mutations_by_tree.sh -> induce_mutations.sh

# mutations.tab: ${MUTFILES}
# 	cat $^ | cut -d' ' -f1,2 | tr ' ' '\t' > $@


# get sites to mutate
#####################

sites_to_mutate.bed: ../e_mel_3/e_mel_3.fa.fai
	bedtools random -l 1 -n 14000 -seed 999 -g $< | bedtools sort -faidx $< > $@

sites_to_mutate.vcf.gz: sites_to_mutate.bed
	$(MAKE) -f randsites.mk $@

sites_to_mutate.txt: sites_to_mutate.vcf.gz
	bcftools query -s GL/M -f '%CHROM\t%POS\t%REF\t%FIRST_ALT\t[%TGT]\n' $< | tr '/' '\t' > $@

#make mut files and split them out according to the tree topology
#################################################################

mut_files/:
	mkdir -p $@

mut_files/mutations.tab: sites_to_mutate.txt mut_files
	cat $< | cut -f 1,2 > $@ 

mut_files/mutations.tab.gz: ../mut_files/mutations.tab mut_files
	bgzip -c <$< >$@

mut_files/mutations.tab.gz.tbi: ../mut_files/mutations.tab.gz mut_files
	tabix -b2 -e2 $<

mut_files/mutfile.txt: sites_to_mutate.txt mut_files/
	cd $(word 2, $^) && cat ../$< | $(SETUPMUTATIONS) && cd ..

$(MUTFILES): mut_files/mutfile.txt mut_files/
	$(SPLITBYTREE) -i $< -o $(word 2, ) -n 1000

#make sed scripts and mutate reads
##################################

sed_scripts/:
	mkdir -p $@

#for all the sed scripts: the r1 and r2 of the sed script is made from: the mutfile with the same stem
$(SEDSCRIPTS): sed_scripts/%.R1.sed sed_scripts/%.R2.sed : $(BAM) -f -r mut_files/% $(REF) sed_scripts/
	$(INDUCEMUTATIONS) -b 

#induce mutations.sh makes 2 files: $(basename $MUTFILE).R1.sed and $(basename $MUTFILE).R2.sed and some log files



#make alignment, but use mutated reads
######################################

alignment.bam: ../e_mel_3/e_mel_3.fa
	${NGMALIGNER} -q '*.fastq.gz' -d ./tmp/ -r $< -o $@ -i ${MUTATEDREADSDIR}

%.dedup.bam: %.bam
	picard MarkDuplicates INPUT=$< OUTPUT=$@ METRICS_FILE=$@.metrics.txt CREATE_INDEX=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

%.dict: %.fa
	picard CreateSequenceDictionary REFERENCE=$< OUTPUT=$@

%.firstcalls.vcf: %.dedup.bam e_mel_3/e_mel_3.fa e_mel_3/e_mel_3.dict
	mkdir -p bqsr; parallel --halt 2 gatk HaplotypeCaller --heterozygosity 0.025 -R $(word 2, $^) -I $< -stand-call-conf 50 -ploidy 2 -L scaffold_{} -O bqsr/{}_$@ ::: $(shell seq 1 11); picard SortVcf $(foreach num, $(shell seq 1 11), I=bqsr/$(num)_$@) O=$@ SEQUENCE_DICTIONARY=$(word 3, $^)

%.recal.table: %.firstcalls.vcf %.dedup.bam e_mel_3/e_mel_3.fa
	gatk BaseRecalibrator -I $(word 2, $^) -R $(word 3, $^) --known-sites $< -O $@

%.recal.bam: %.recal.table %.dedup.bam e_mel_3/e_mel_3.fa
	gatk ApplyBQSR -I $(word 2, $^) -R $(word 3, $^) --bqsr-recal-file $< -O $@

