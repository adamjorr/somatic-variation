#!usr/bin/bash

#bwa_aligner.sh reference.fasta
#Takes reference and aligns all .fastq files in any subdirectory.
#This is a modified version of align_fastq which uses bwa instead of bowtie.
USAGE="Usage: $0 [-t THREADS] [-p RG_PLATFORM] -r reference.fa -o out.bam data/"




#Here are some things you might want to change:
RGPL=ILLUMINA #We assume Illumina; if we're wrong, change it here.
CORES=48
REFERENCEFILE=""
OUTNAME=""

while getopts :t:p:r:o:h opt; do
	case $opt in
		t)
			CORES=$OPTARG
			;;
		p)
			RGPL=$OPTARG
			;;
		r)
			REFERENCEFILE=$OPTARG
			;;
		o)
			OUTNAME=$OPTARG
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

shift $((OPTIND-1)) # get operands
if [ $# -ne 1 ]; then			#if we forget arguments
	echo $USAGE	#remind us
	exit 1				#and exit with error
fi

if [ $REFERENCEFILE == "" ]; then
	echo $USAGE
	echo "Reference file required."
	exit 1
fi

if [ $OUTNAME == "" ]; then
	echo $USAGE
	echo "Output file name required."
	exit 1
fi

#Some variables
DATADIR=$1
FASTQFILES=$(find $DATADIR -name '*R1*.fastq') || exit
NUMFILES=( $FASTQFILES )
NUMFILES=${#NUMFILES[@]}

if [[ ! -e ${REFERENCEFILE}.amb || ! -e ${REFERENCEFILE}.ann || ! -e ${REFERENCEFILE}.bwt || ! -e ${REFERENCEFILE}.pac || ! -e ${REFERENCEFILE}.sa ]]; then
	bwa index ${REFERENCEFILE}
fi

echo Making BAM files . . . >&2
#Make bamfiles from the FASTQs
PROGRESS=0
for F in $FASTQFILES; do
	((PROGRESS++))
	echo Progress: $PROGRESS / $NUMFILES >&2
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

echo Merging . . . >&2
samtools merge -@ $CORES $OUTNAME $BAMS || exit

#Now clean up
rm $BAMS || exit

exit 0
