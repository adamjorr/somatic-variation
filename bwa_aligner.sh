#!usr/bin/bash

#bwa_aligner.sh reference.fasta
#Takes reference and aligns all .fastq files in any subdirectory.
#This is a modified version of align_fastq which uses bwa instead of bowtie.
USAGE="Usage: $0 [-t THREADS] [-p RG_PLATFORM] reference.fa data/ out.bam"


if [ $# -ne 3 ]; then			#if we forget arguments
	echo $USAGE	#remind us
	exit 1				#and exit with error
fi

#Here are some things you might want to change:
RGPL=ILLUMINA #We assume Illumina; if we're wrong, change it here.
CORES=48

while getopts :t:p:h opt; do
	shift $((OPTIND-1))
	case $opt in
		t)
			CORES=$OPTARG
			;;
		p)
			RGPL=$OPTARG
			;;
		h)
			echo $USAGE >&2
			exit 1
			;;
	    \?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done



#Some variables
REFERENCEFILE=$1
DATADIR=$2
OUTNAME=$3
FASTQFILES=$(find $DATADIR -name '*R1*.fastq') || exit

if [[ ! -e ${REFERENCEFILE}.amb || ! -e ${REFERENCEFILE}.ann || ! -e ${REFERENCEFILE}.bwt || ! -e ${REFERENCEFILE}.pac || ! -e ${REFERENCEFILE}.sa ]]; then
	bwa index ${REFERENCEFILE}
fi

echo Making BAM files . . .
#Make bamfiles from the FASTQs
for F in $FASTQFILES; do
	BASEFNAME=$(basename $F) || exit
	BAMS=$(echo $BAMS ${BASEFNAME%R1*}.bam) || exit
	if [ -e ${BASEFNAME%R1*}.bam ]; then samtools quickcheck ${BASEFNAME%R1*}.bam || rm ${BASEFNAME%R1*}.bam; fi
	if [ ! -e ${BASEFNAME%R1*}.bam ]; then
		RGPU=$(head -n 1 $F | cut -d: -f3,4 --output-delimiter=.) || exit
		RGLB=$(expr $F : '.*\(M[0-9]*[abc]\)') || exit
		RGSM=$(expr $F : '.*\(M[0-9]*[abc]\)') || exit
		bwa mem -t ${CORES} -M -R '@RG\tID:'${RGSM}'\tPL:'${RGPL}'\tPU:'${RGPU}'\tLB:'${RGLB}'\tSM:'${RGSM} $REFERENCEFILE $F ${F/R1/R2} > ${BASEFNAME%R1*}.sam
		samtools sort -@ $CORES -o ${BASEFNAME%R1*}.bam -m 2G -T tmp ${BASEFNAME%R1*}.sam || exit
		rm ${BASEFNAME%R1*}.sam || exit
	fi
done

echo Merging . . .
samtools merge -@ $CORES $OUTNAME $BAMS || exit

#Now clean up
rm $BAMS || exit

exit 0
