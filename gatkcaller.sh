#!usr/bin/bash
##TODO: REMOVE ALIGNMENT COMPONENTS, UPDATE DOCS + VARIABLES
#aligntoref.sh reference.fasta data.bam
#Takes reference and aligns all .fastq files in any subdirectory and calls SNPs with GATK.
#requires samtools, picard-tools, the GATK, and bowtie2

#http://tldp.org/LDP/abs/html/string-manipulation.html is a great guide for manipulating strings in bash

if [ $# -ne 1 ]; then			#if we forget arguments
	echo "Usage: $0 reference.fasta data.bam"	#remind us
	exit 1				#and exit with error
fi

#Here are some things you might want to change:
PICARD="picard" #How do I call picard on this system?
GATK=~/bin/GenomeAnalysisTK.jar #Location of your GATK jar

#Some variables
REFERENCEFILE=$1
FILEIN=$2
CORES=16
HALFCORES=$((CORES / 2))

###Below uses GATK to do some analysis.



$PICARD MarkDuplicates INPUT=$FILEIN OUTPUT=dedup.bam METRICS_FILE=metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 || exit
samtools faidx $REFERENCEFILE || exit
$PICARD CreateSequenceDictionary REFERENCE=${REFERENCEFILE} OUTPUT=${REFERENCEFILE%.*}.dict
$PICARD BuildBamIndex INPUT=dedup.bam

#Now we use GATK to recalibrate our quality scores and give us a VCF.
java -jar ${GATK} -T RealignerTargetCreator -nt $CORES -R ${REFERENCEFILE} -I dedup.bam -o forIndelAligner.intervals || exit
java -jar ${GATK} -T IndelRealigner -R ${REFERENCEFILE} -I dedup.bam -targetIntervals forIndelAligner.intervals -o realigned.bam || exit
java -jar ${GATK} -T UnifiedGenotyper -I realigned.bam -R ${REFERENCEFILE} -stand_call_conf 50 -stand_emit_conf 50 -ploidy 2 -glm BOTH -o first-calls.vcf || exit
java -jar ${GATK} -T BaseRecalibrator -nct $CORES -I realigned.bam -R ${REFERENCEFILE} --knownSites first-calls.vcf -o recal_data.table || exit
java -jar ${GATK} -T PrintReads -nct $CORES -I realigned.bam -R ${REFERENCEFILE} -BQSR recal_data.table -EOQ -o recal.bam || exit
java -jar ${GATK} -T UnifiedGenotyper -nt $HALFCORES -nct $HALFCORES -I recal.bam -R ${REFERENCEFILE} -ploidy 2 -glm BOTH -o var-calls.vcf || exit

#Clean up some things
rm ${REFERENCEFILE}.fai ${REFERENCEFILE%.*}.dict metrics.txt aligner* forIndelAligner.intervals realigned.bam first-calls.vcf recal_data.table dedup.bam

exit 0
