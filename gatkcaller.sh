#!usr/bin/bash

#aligntoref.sh reference.fasta
#Takes reference and aligns all .fastq files in any subdirectory and calls SNPs with GATK.
#requires samtools, picard-tools, the GATK, and bowtie2

#http://tldp.org/LDP/abs/html/string-manipulation.html is a great guide for manipulating strings in bash

if [ $# -ne 1 ]; then			#if we forget arguments
	echo "Usage: $0 reference.fasta"	#remind us
	exit 1				#and exit with error
fi

#Here are some things you might want to change:
PLATFORM=ILLUMINA #We assume Illumina; if we're wrong, change it here.
PICARD="picard" #How do I call picard on this system?
GATK=~/bin/GenomeAnalysisTK.jar #Location of your GATK jar

#Some variables
REFERENCEFILE=$1
FASTQFILES=$(find ./data/ -name '*R1*.fastq')

#Build fasta index
if [ ! -e aligner.1.bt2 ] || [ ! -e aligner.rev.2.bt2 ]; then
	bowtie2-build $REFERENCEFILE aligner >/dev/null || exit
fi

#Make bamfiles from the FASTQs
for F in $FASTQFILES; do
	BASEFNAME=$(basename $F)
	if [ ! -e ${BASEFNAME%R1*}.bam ]; then
		RGPU=$(head -n 1 $F | cut -d: -f3,4 --output-delimiter=.)
		RGLB=$(expr $F : '\(M[0-9]*\)')
		RGSM=$(expr $F : '\(M[0-9]*[abc]\)')
		bowtie2 -x aligner --phred33 --rg-id ${RGSM} --rg PL:${PLATFORM} --rg PU:${RGPU} --rg LB:${RGLB} --rg SM:${RGSM} -1 $F -2 ${F/R1/R2} -S ${BASEFNAME%R1*}.sam >/dev/null || exit
		samtools view -bh -o ${BASEFNAME%R1*}.bam ${BASEFNAME%R1*}.sam || exit
		rm ${BASEFNAME%R1*}.sam || exit
	fi
done

for F in $(find . -name '*.bam'); do
	INPUTS=$(echo $INPUTS INPUT=$F)
done

#Now we merge the files, mark duplicates, make a sequence dictionary and a bam index.
$PICARD MergeSamFiles $INPUTS OUTPUT=data.bam USE_THREADING=true || exit
$PICARD MarkDuplicates INPUT=data.bam OUTPUT=dedup.bam METRICS_FILE=metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 || exit
samtools faidx $REFERENCEFILE || exit
$PICARD CreateSequenceDictionary REFERENCE=${REFERENCEFILE} OUTPUT=${REFERENCEFILE%.*}.dict
$PICARD BuildBamIndex INPUT=dedup.bam

#Now we use GATK to recalibrate our quality scores and give us a VCF.
java -jar ${GATK} -T RealignerTargetCreator -R ${REFERENCEFILE} -I dedup.bam -o forIndelAligner.intervals || exit
java -jar ${GATK} -T IndelRealigner -R ${REFERENCEFILE} -I dedup.bam -targetIntervals forIndelAligner.intervals -o realigned.bam || exit
java -jar ${GATK} -T UnifiedGenotyper -I realigned.bam -R ${REFERENCEFILE} -stand_call_conf 50 -stand_emit_conf 50 -ploidy 2 -glm BOTH -o first-calls.vcf || exit
java -jar ${GATK} -T BaseRecalibrator -I realigned.bam -R ${REFERENCEFILE} --knownSites first-calls.vcf -o recal_data.table || exit
java -jar ${GATK} -T PrintReads -I realigned.bam -R ${REFERENCEFILE} -BQSR recal_data.table -EOQ -o recal.bam || exit
java -jar ${GATK} -T UnifiedGenotyper -I recal.bam -R ${REFERENCEFILE} -ploidy 2 -glm BOTH -o var-calls.vcf || exit

#Clean up some things
rm ${REFERENCEFILE}.fai ${REFERENCEFILE%.*}.dict metrics.txt aligner* forIndelAligner.intervals realigned.bam first-calls.vcf recal_data.table

exit
