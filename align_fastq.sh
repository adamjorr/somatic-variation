#!usr/bin/bash

#aligntoref.sh reference.fasta
#Takes reference and aligns all .fastq files in any subdirectory.
#requires samtools, picard-tools, and bowtie2

#http://tldp.org/LDP/abs/html/string-manipulation.html is a great guide for manipulating strings in bash

if [ $# -ne 1 ]; then			#if we forget arguments
	echo "Usage: $0 reference.fasta"	#remind us
	exit 1				#and exit with error
fi

#Here are some things you might want to change:
PLATFORM=ILLUMINA #We assume Illumina; if we're wrong, change it here.
PICARD="picard" #How do I call picard on this system?
CORES=16

#Some variables
REFERENCEFILE=$1
FASTQFILES=$(find ./data/ -name '*R1*.fastq') || exit

#Build fasta index
if [ ! -e aligner.1.bt2 ] || [ ! -e aligner.rev.2.bt2 ]; then
	bowtie2-build $REFERENCEFILE aligner || exit
fi

echo Making BAM files . . .
#Make bamfiles from the FASTQs
for F in $FASTQFILES; do
	BASEFNAME=$(basename $F) || exit
	BAMS=$(echo $BAMS ${BASEFNAME%R1*}.bam) || exit
	if [ ! -e ${BASEFNAME%R1*}.bam ]; then
		RGPU=$(head -n 1 $F | cut -d: -f3,4 --output-delimiter=.) || exit
		RGLB=$(expr $F : '\(M[0-9]*[abc]\)') || exit
		RGSM=$(expr $F : '\(M[0-9]*[abc]\)') || exit
		bowtie2 -p $CORES -x aligner --phred33 --rg-id ${RGSM} --rg PL:${PLATFORM} --rg PU:${RGPU} --rg LB:${RGLB} --rg SM:${RGSM} -1 $F -2 ${F/R1/R2} -S ${BASEFNAME%R1*}.sam || exit
		samtools sort -@ $CORES -o ${BASEFNAME%R1*}.bam -n -T tmp ${BASEFNAME%R1*}.sam || exit
		rm ${BASEFNAME%R1*}.sam || exit
	fi
done

echo Merging . . .
#Now we merge the files
#$PICARD MergeSamFiles $INPUTS OUTPUT=data.bam USE_THREADING=true || exit
samtools merge -@ $CORES -n data.bam $BAMS || exit

#Now clean up
rm $BAMS aligner* || exit

exit 0