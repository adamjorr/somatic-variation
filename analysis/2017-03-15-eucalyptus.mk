READSFOLDER=./data/
CLEANFOLDER=./cleaned_reads/sliced/
BWALIGNER=~/bin/zypy/bwa_aligner.sh
NGMALIGNER=~/bin/zypy/ngm_aligner.sh
CREATECONSENSUS=~/bin/zypy/create_consensus.sh
GATKCALLER=~/bin/zypy/gatkcaller.sh
READS=$(shell find ${READSFOLDER} -name "*.fastq" -and -name "*R1*")
READNAMES=$(notdir $(READS))
CORREADS=$(subst .fastq,.cor.fq,$(READNAMES))
SLICEDNAMES=$(subst .cor.fq,.cor_sliced.fastq,$(CORREADS))
SLICEDREADS=$(addprefix $(CLEANFOLDER),$(SLICEDNAMES))
CLEANREADS=./clean_reads.sh
REF=./ref.fa


all: filtered_bed_excluded.nwk

filtered_bed_included.fa: filtered.vcf
	cat $< | filt_with_replicates.pl -s -g 3 | vcf2fa.sh | diploidify.py -v > $@

filtered_bed_excluded.nwk: filtered_bed_excluded.vcf
	vcf2tree.sh -t 16 -i $< -o $@ -g 3

filtered_bed_excluded.fa: filtered_bed_excluded.vcf
	 cat $< | filt_with_replicates.pl -s -g 3 | vcf2fa.sh | diploidify.py -v > $@

filtered_bed_excluded_no_s.nwk: filtered_bed_excluded_no_s.fa
	raxml -T 8 -f a -s $< -n nwk -m ASC_GTRGAMMA -w /tmp/ --asc-corr=lewis -p 12345 -x 12345 -# 100 && cat /tmp/RAxML_bestTree.nwk >$@

filtered_bed_excluded_no_s.fa: filtered_bed_excluded.vcf
	cat $< | filt_with_replicates.pl -g 3 | vcf2fa.sh | diploidify.py -v > $@

filtered_bed_excluded_majority.nwk: filtered_bed_excluded_majority.fa
	raxml -T 16 -f a -s $< -n nwk -m ASC_GTRGAMMA -w /tmp/ --asc-corr=lewis -p 12345 -x 12345 -# 100 && cat /tmp/RAxML_bestTree.nwk >$@

filtered_bed_excluded_majority.fa: filtered_bed_excluded.vcf
	cat $< | filt_with_replicates.pl --majority -s -g 3 | vcf2fa.sh | diploidify.py -v > $@

dng_calls.vcf: dng/sampleM.ped filtered_bed_excluded.ad
	dng call --ped=$< --mu=0 --mu-library=0 --mu-somatic=1e-7 --theta=0.026 --lib-error=9.319e-04 --lib-overdisp-hom=0.007 --lib-overdisp-het=0.306 --lib-bias=1 --min-prob=0 --ref-bias-hom=0.957 --ref-bias-het=0.856 --ref-bias-hap=0 --input=$(word 2, $^) --output=$@

filtered_bed_excluded.ad: filtered_bed_excluded.sites alignment.bam
	dng pileup -o $@ -r @$< -f e_mel_3/e_mel_3.fa $(word 2, $^)

filtered_bed_excluded.sites: filtered_bed_excluded.vcf
	bcftools view -H $< | awk '{print $$1":"$$2'} > $@

filtered_bed_excluded.vcf: filtered.vcf liftover/e_mel_3_repeatmask.bed
	vcftools --vcf $< --exclude-bed $(word 2, $^) --remove-filtered-all --recode --recode-INFO-all --stdout >$@

filtered.vcf: var-calls-first11.vcf
	bcftools filter -g 50 -i 'DP <= 500 && ExcessHet <=40' $< | bcftools view -m2 -M2 -v snps -o $@

var-calls-first11.vcf: var-calls.vcf
	gzip -c $< | bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8,scaffold_9,scaffold_10,scaffold_11 -o $@

tree.nwk: var-calls.vcf
	vcf2tree.sh -t 16 -i $< -o $@ -g 3

var-calls.vcf: e_mel_3/e_mel_3.fa alignment.bam
	${GATKCALLER} -d ./tmp/ -o $@ -r $< -i $(word 2, $^)

alignment.bam: e_mel_3/e_mel_3.fa
	${NGMALIGNER} -d ./tmp/ -r $< -o $@ -i ${READSFOLDER}

.PRECIOUS: alignment.bam %.dedup.bam %.realigned.bam %.realigned.interval_list %.scattered.intervals/ %.firstcalls.vcf %.recal.table %.recal.bam

%.dedup.bam: %.bam
	picard MarkDuplicates INPUT=$< OUTPUT=$@ METRICS_FILE=$@.metrics.txt CREATE_INDEX=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

%.realigned.interval_list: %.dedup.bam e_mel_3/e_mel_3.fa
	java -jar ~/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 48 -R $(word 2, $^) -I $< -o $@

%.scattered.intervals/: %.realigned.interval_list
	mkdir -p $@; picard IntervalListTools I=$< SCATTER_COUNT=32 O=$@

%.realigned.bam: %.dedup.bam e_mel_3/e_mel_3.fa %.realigned.interval_list
	mkdir -p realn; parallel --halt 2 java -Xmx100G -jar ~/bin/GenomeAnalysisTK.jar -T IndelRealigner -R $(word 2, $^) -I $< -targetIntervals $(word 3, $^) -L scaffold_{} -o realn/{}_$@ ::: $(shell seq 1 11); picard GatherBamFiles $(foreach num, $(shell seq 1 11), I=realn/$(num)_$@) O=$@

%.firstcalls.vcf: %.realigned.bam e_mel_3/e_mel_3.fa e_mel_3/e_mel_3.dict
	mkdir -p bqsr; parallel --halt 2 gatk HaplotypeCaller --heterozygosity 0.025 -R $(word 2, $^) -I $< -stand-call-conf 50 -ploidy 2 -L scaffold_{} -O bqsr/{}_$@ ::: $(shell seq 1 11); picard SortVcf $(foreach num, $(shell seq 1 11), I=bqsr/$(num)_$@) O=$@ SEQUENCE_DICTIONARY=$(word 3, $^)

%.recal.table: %.firstcalls.vcf %.realigned.bam e_mel_3/e_mel_3.fa
	gatk BaseRecalibrator -I $(word 2, $^) -R $(word 3, $^) --known-sites $< -O $@

%.recal.bam: %.recal.table %.realigned.bam e_mel_3/e_mel_3.fa
	gatk ApplyBQSR -I $(word 2, $^) -R $(word 3, $^) --bqsr-recal-file $< -O $@

liftover/e_mel_1_repeatmask.gff3: liftover/egrandis_201_repeatmask.gff3 e_mel_1/e_mel_1.fa
	liftOver -gff $< e_mel_1.fa.chain $@ /dev/null

liftover/e_mel_2_repeatmask.gff3: liftover/e_mel_1_repeatmask.gff3 e_mel_2/e_mel_2.fa
	liftOver -gff $< e_mel_2.fa.chain $@ /dev/null

liftover/e_mel_3_repeatmask.gff3: liftover/e_mel_2_repeatmask.gff3 e_mel_3/e_mel_3.fa
	liftOver -gff $< e_mel_3.fa.chain $@ /dev/null

liftover/e_mel_3_repeatmask.bed: liftover/e_mel_3_repeatmask.gff3
	gff2bed <$< >$@

liftover/e_mel_1_genes.gff3: liftover/egrandis_201_genes.gff3 e_mel_1/e_mel_1.fa
	liftOver -gff $< e_mel_1.fa.chain $@ /dev/null

liftover/e_mel_2_genes.gff3: liftover/e_mel_1_genes.gff3 e_mel_2/e_mel_2.fa
	liftOver -gff $< e_mel_2.fa.chain $@ /dev/null

liftover/e_mel_3_genes.gff3: liftover/e_mel_2_genes.gff3 e_mel_3/e_mel_3.fa
	liftOver -gff $< e_mel_3.fa.chain $@ /dev/null

overlapping_randsites.bam: e_mel_3_randsites.bed alignment.bam
	samtools view -h -b -o $@ -L $< $(word 2, $^)

e_mel_3_randsites.bed: e_mel_3/e_mel_3.fa.fai
	bedtools random -l 1 -n 150000 -seed 999 -g $< | bedtools sort -faidx $< > $@

e_mel_3_randsites_norepeats.bed: e_mel_3_randsites.bed #liftover/e_mel_3_repeatmask.bed
	bedtools intersect -v -a e_mel_3_randsites.bed -b liftover/e_mel_3_repeatmask.bed > $@

e_mel_3_randsites.regions: e_mel_3_randsites.bed
	awk '{print $$1":"$$3'} < $< > $@

e_mel_3_randsites_norepeats.regions: e_mel_3_randsites_norepeats.bed
	awk '{print $$1":"$$3'} < $< > $@

overlapping_randsites.ad: e_mel_3_randsites.regions
	dng pileup -Q 13 -q 3 -o $@ -r @e_mel_3_randsites.regions -f e_mel_3/e_mel_3.fa overlapping_randsites.bam

overlapping_randsites_norepeats.ad: e_mel_3_randsites_norepeats.regions
	dng pileup -Q 13 -q 3 -o $@ -r @e_mel_3_randsites_norepeats.regions -f e_mel_3/e_mel_3.fa overlapping_randsites.bam

e_mel_3/e_mel_3.fa: e_mel_3/e_mel_3.bam e_mel_2/e_mel_2.fa
	${CREATECONSENSUS} -o $@ -r $(word 2, $^) -i $<

e_mel_2/e_mel_2.fa: e_mel_2/e_mel_2.bam e_mel_1/e_mel_1.fa
	${CREATECONSENSUS} -o $@ -r $(word 2, $^) -i $<

e_mel_1/e_mel_1.fa: e_mel_1/e_mel_1.bam ${REF}
	${CREATECONSENSUS} -o $@ -r $(word 2, $^) -i $<

e_mel_3/e_mel_3.bam: e_mel_2/e_mel_2.fa $(SLICEDREADS)
	mkdir -p e_mel_3
	${NGMALIGNER} -s 0.3 -d ./tmp/ -r $< -o $@ -i ${CLEANFOLDER}

e_mel_2/e_mel_2.bam: e_mel_1/e_mel_1.fa $(SLICEDREADS)
	mkdir -p e_mel_2
	${NGMALIGNER} -s 0.3 -d ./tmp/ -r $< -o $@ -i ${CLEANFOLDER}

e_mel_1/e_mel_1.bam: ${REF} $(SLICEDREADS)
	mkdir -p e_mel_1
	${NGMALIGNER} -s 0.3 -d ./tmp/ -r $< -o $@ -i ${CLEANFOLDER}

$(SLICEDREADS) : $(READS)
	${CLEANREADS} -i ${READSFOLDER}


e_mel_3_cds.bed: liftover/e_mel_3_genes.gff3
	cat $< | awk '$$3 == "CDS" {print}' | gff2bed > $@

degenerate/e_mel_3_degenerate.sites: e_mel_3_cds.bed e_mel_3/e_mel_3.fa
	python2 ~/bin/zypy/misc/degenerate.py -i $< -r $(word 2, $^) -o $@

degenerate/e_mel_3_degenerate.bed: degenerate/e_mel_3_degenerate.sites
	cat $< | tr : '\t' > $@

# Reed's Rules
e_mel_3/e_mel_3.bed: e_mel_3/e_mel_3.fa.fai
	awk 'BEGIN {FS="\t"}; {print $$1 FS "0" FS $$2}' $< > $@

e_mel_3/e_mel_3_norepeats.bed: e_mel_3/e_mel_3.bed
	bedtools subtract -a $< -b liftover/e_mel_3_repeatmask.bed \
	| bedtools merge -d 10 \
	| awk '$$3-$$2 >= 200' \
	> $@

e_mel_3/e_mel_3_norepeats.regions: e_mel_3/e_mel_3_norepeats.bed
	awk '{print $$1 ":" $$2+1 "-" $$3'} < $< > $@

alignment_norepeats.ad:
	dng pileup -Q 13 -q 3 --max-dp 500 --min-dp 16 -o $@ -r @e_mel_3_norepeats.regions -f e_mel_3/e_mel_3.fa alignment.bam
