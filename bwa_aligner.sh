#!usr/bin/bash

#bwa_aligner.sh reference.fasta
#Takes reference and aligns all .fastq files in any subdirectory.
#This is a modified version of align_fastq which uses bwa instead of bowtie.

if [ $# -ne 3 ]; then			#if we forget arguments
	echo "Usage: $0 reference.fa data/ out.bam"	#remind us
	exit 1				#and exit with error
fi

#Here are some things you might want to change:
PLATFORM=ILLUMINA #We assume Illumina; if we're wrong, change it here.
CORES=16

#Some variables
REFERENCEFILE=$1
DATADIR=$2
OUTNAME=$3
FASTQFILES=$(find $DATADIR -name '*R1*.fastq') || exit

echo Making BAM files . . .
#Make bamfiles from the FASTQs
for F in $FASTQFILES; do
	BASEFNAME=$(basename $F) || exit
	BAMS=$(echo $BAMS ${BASEFNAME%R1*}.bam) || exit
	if [ ! -e ${BASEFNAME%R1*}.bam ]; then
		RGPU=$(head -n 1 $F | cut -d: -f3,4 --output-delimiter=.) || exit
		RGLB=$(expr $F : '.*\(M[0-9]*[abc]\)') || exit
		RGSM=$(expr $F : '.*\(M[0-9]*[abc]\)') || exit
		bwa mem -t 16 -M -R '@RG\tID:'${RGSM}'\tPL:'${RGPL}'\tPU:'${RGPU}'\tLB:'${RGLB}'\tSM:'${RGSM} $REFERENCEFILE $F ${F/R1/R2} > ${BASEFNAME%R1*}.sam
		samtools sort -@ $CORES -o ${BASEFNAME%R1*}.bam -n -T tmp ${BASEFNAME%R1*}.sam || exit
		rm ${BASEFNAME%R1*}.sam || exit
	fi
done

echo Merging . . .
samtools merge -@ $CORES -n $OUTNAME $BAMS || exit

#Now clean up
rm $BAMS || exit

exit 0
