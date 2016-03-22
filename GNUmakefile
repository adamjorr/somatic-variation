NUMBERS := $(shell seq 2 $(MAKECMDGOALS))
NUMFILES := $(addsuffix .bam,$(NUMBERS))
SCRIPTDIR :=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

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
	echo making 1.bam
	echo bash $(SCRIPTDIR)/bwa_aligner.sh $(REFFASTA) $(DATADIR) 1_bwa.bam
	samtools view -@ 16 -h -b -U 1_bwa_unmapped.bam -o 1_bwa_mapped.bam -@ 16 -f 2 1_bwa.bam
	rm 1_bwa.bam
	echo bash $(SCRIPTDIR)/stampy_realigner.sh $(REFFASTA) 1_bwa_unmapped.bam 1_stampy.bam

$(NUMBERS) : % : %.bam

.SECONDEXPANSION:

$(NUMFILES) : %.bam : $$(shell expr $$* - 1).bam
	echo making $@ from $<

