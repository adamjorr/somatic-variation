NUMBERS := $(shell seq 2 $(MAKECMDGOALS))
NUMFILES := $(addsuffix .bam,$(NUMBERS))
NUMREFS := $(addsuffix .fa,$(NUMBERS))
SCRIPTDIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

QUALITY := 13
CORES := 16
GATK := java -jar ~/bin/GenomeAnalysisTK.jar
SEQTK := ~/bin/seqtk/seqtk

ifeq ($(strip $(MAKECMDGOALS)),0)
$(error Call make with an argument > 0)
else ifeq ($(strip $(MAKECMDGOALS)),)
$(error Call make with an argument > 0)
endif
ifeq ($(strip $(REFFASTA)),)
$(error Call make with REFFASTA=ref.fa)
endif
ifeq ($(strip $(DATADIR)),)
$(error Call make with DATADIR=../data/)
endif


1 : 1.bam

1.bam :
	bash $(SCRIPTDIR)/bwa_aligner.sh $(REFFASTA) $(DATADIR) 1_bwa.bam
	samtools view -@ $(CORES) -b -h -q $(QUALITY) -f 2 -o 1_bwa_mapped.bam -U 1_bwa_unmapped.bam 1_bwa.bam
	rm 1_bwa.bam
	bash $(SCRIPTDIR)/stampy_realigner.sh $(REFFASTA) 1_bwa_unmapped.bam 1_stampy.bam
	rm 1_bwa_unmapped.bam
	samtools merge -@ $(CORES) -n -c -p 1_merged.bam 1_bwa_mapped.bam 1_stampy.bam
	rm 1_bwa_mapped.bam 1_stampy.bam
	samtools sort -@ $(CORES) -m 2G -o 1_sorted.bam -O bam -T tmp 1_merged.bam
	rm 1_merged.bam
	$(GATK) -T RealignerTargetCreator -R $(REFFASTA) -I 1_merged.bam -o realignment_targets.list
	$(GATK) -T IndelRealigner -R $(REFFASTA) -I 1_merged.bam -targetIntervals realignment_targets.list -o 1.bam
	rm realignment_targets.list 1_merged.bam
	samtools index -@ $(CORES) 1.bam

1.fa : 1.bam
	samtools mpileup -uf $(REFFASTA) 1.bam | bcftools call -c - | vcfutils.pl vcf2fq | $(SEQTK) -a - > 1.fa

$(NUMBERS) : % : %.bam

.SECONDEXPANSION:

$(NUMFILES) : %.bam : $$(shell expr $$* - 1).fa
	bash $(SCRIPTDIR)/bwa_aligner.sh $< $(DATADIR) $*_bwa.bam
	samtools view -@ $(CORES) -b -h -q $(QUALITY) -f 2 -o $*_bwa_mapped.bam -U $*_bwa_unmapped.bam $*_bwa.bam
	rm $*_bwa.bam
	bash $(SCRIPTDIR)/stampy_realigner.sh $(REFFASTA) $*_bwa_unmapped.bam $*_stampy.bam
	rm $*_bwa_unmapped.bam
	samtools merge -@ $(CORES) -n -c -p $*_merged.bam $*_bwa_mapped.bam $*_stampy.bam
	rm $*_bwa_mapped.bam $*_stampy.bam
	samtools sort -@ $(CORES) -m 2G -o $*_sorted.bam -O bam -T tmp $*_merged.bam
	rm $*_merged.bam
	$(GATK) -T RealignerTargetCreator -R $< -I $*_merged.bam -o realignment_targets.list
	$(GATK) -T IndelRealigner -R $(REFFASTA) -I $*_merged.bam -targetIntervals realignment_targets.list -o $@
	rm realignment_targets.list $*_merged.bam
	samtools index -@ $(CORES) $@

$(NUMREFS) : %.fa : $$(shell expr $$* - 1).fa %.bam
	samtools mpileup -uf $< $*.bam | bcftools call -c - | vcfutils.pl vcf2fq | $(SEQTK) -a - > $@
