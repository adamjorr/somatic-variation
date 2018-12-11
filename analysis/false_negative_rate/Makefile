SCRIPTDIR=../scripts
NGMALIGNER=$(SCRIPTDIR)/ngm_aligner.sh

READSFOLDER=24_1_2018-mutated_fastq
MUTFILES=$(shell find 24_1_2018-mutated_fastq -name "mut_*.txt")

#make alignment, but use mutated reads

alignment.bam: ../e_mel_3/e_mel_3.fa
	${NGMALIGNER} -q '*.fastq.gz' -d ./tmp/ -r $< -o $@ -i ${READSFOLDER}

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


mutations.tab: ${MUTFILES}
	cat $^ | cut -d' ' -f1,2 | tr ' ' '\t' > $@

mutfile.txt: ${MUTFILES}
	cat $^ | sort | uniq > $@

