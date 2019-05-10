READSFOLDER=./data/
REF=./ref.fa
REPEATMASK=./mask.gff3
CLEANFOLDER=./cleaned_reads/sliced/
NGMALIGNER=ngm_aligner.sh
CREATECONSENSUS=create_consensus.sh
GATKCALLER=gatkcaller.sh
CLEANREADS=clean_reads.sh
READS=$(shell find ${READSFOLDER} -name "*.fastq" -and -name "*R1*")
READNAMES=$(notdir $(READS))
CORREADS=$(subst .fastq,.cor.fq,$(READNAMES))
SLICEDNAMES=$(subst .cor.fq,.cor_sliced.fastq,$(CORREADS))
SLICEDREADS=$(addprefix $(CLEANFOLDER),$(SLICEDNAMES))

all: filtered_bed_excluded.nwk

$(SLICEDREADS) : $(READS)
    ${CLEANREADS} -i ${READSFOLDER}

e_mel_1/e_mel_1.bam: ${REF} $(SLICEDREADS)
    mkdir -p e_mel_1
    ${NGMALIGNER} -s 0.3 -d ./tmp/ -r $< -o $@ -i ${CLEANFOLDER}

e_mel_1/e_mel_1.fa: e_mel_1/e_mel_1.bam ${REF}
    ${CREATECONSENSUS} -o $@ -r $(word 2, $^) -i $<

liftover/e_mel_1_repeatmask.gff3: $(REPEATMASK) e_mel_1/e_mel_1.fa
    liftOver -gff $< e_mel_1.fa.chain $@ /dev/null

e_mel_2/e_mel_2.bam: e_mel_1/e_mel_1.fa $(SLICEDREADS)
    mkdir -p e_mel_2
    ${NGMALIGNER} -s 0.3 -d ./tmp/ -r $< -o $@ -i ${CLEANFOLDER}

e_mel_2/e_mel_2.fa: e_mel_2/e_mel_2.bam e_mel_1/e_mel_1.fa
    ${CREATECONSENSUS} -o $@ -r $(word 2, $^) -i $<

liftover/e_mel_2_repeatmask.gff3: liftover/e_mel_1_repeatmask.gff3 e_mel_2/e_mel_2.fa
    liftOver -gff $< e_mel_2.fa.chain $@ /dev/null

e_mel_3/e_mel_3.bam: e_mel_2/e_mel_2.fa $(SLICEDREADS)
    mkdir -p e_mel_3
    ${NGMALIGNER} -s 0.3 -d ./tmp/ -r $< -o $@ -i ${CLEANFOLDER}

e_mel_3/e_mel_3.fa: e_mel_3/e_mel_3.bam e_mel_2/e_mel_2.fa
    ${CREATECONSENSUS} -o $@ -r $(word 2, $^) -i $<

liftover/e_mel_3_repeatmask.gff3: liftover/e_mel_2_repeatmask.gff3 e_mel_3/e_mel_3.fa
    liftOver -gff $< e_mel_3.fa.chain $@ /dev/null

liftover/e_mel_3_repeatmask.bed: liftover/e_mel_3_repeatmask.gff3
    gff2bed <$< >$@

alignment.bam: e_mel_3/e_mel_3.fa
    ${NGMALIGNER} -d ./tmp/ -r $< -o $@ -i ${READSFOLDER}

var-calls.vcf: e_mel_3/e_mel_3.fa alignment.bam
    ${GATKCALLER} -d ./tmp/ -o $@ -r $< -i $(word 2, $^)

var-calls.vcf.gz: var-calls.vcf
    bgzip var-calls.vcf

var-calls.vcf.gz.csi: var-calls.vcf.gz
    bcftools index $<

var-calls-first11.vcf: var-calls.vcf.gz
    bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8,scaffold_9,scaffold_10,scaffold_11 -o $@ $<

filtered.vcf: var-calls-first11.vcf
    bcftools filter -g 50 -i 'DP <= 500 && ExcessHet <=40' $< | bcftools view -m2 -M2 -v snps -o $@

filtered_bed_excluded.vcf: filtered.vcf liftover/e_mel_3_repeatmask.bed
    vcftools --vcf $< --exclude-bed $(word 2, $^) --remove-filtered-all --recode --recode-INFO-all --stdout >$@

filtered_bed_excluded.nwk: filtered_bed_excluded.vcf
    vcf2tree.sh -t 16 -i $< -o $@ -g 3

