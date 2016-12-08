#!usr/bin/bash

#bwa_aligner.sh reference.fasta
#Takes reference and aligns all .fastq files in any subdirectory.
#This is a modified version of align_fastq which uses bwa instead of bowtie.
USAGE="Usage: $0 [-t THREADS] [-p RG_PLATFORM] [-q FILEPATTERN] [-1 FIRSTMATE] [-2 SECONDMATE] -r reference.fa -o out.bam data/"




#Here are some things you might want to change:
RGPL=ILLUMINA #We assume Illumina; if we're wrong, change it here.
CORES=48
FILE_PATTERN='*.fastq' #pattern to find first set of reads
SEARCH_STRING='R1' #pattern for search/replace to find second set of reads
REPLACE_STRING='R2' #pattern for search/replace to substitute second set of reads
REFERENCEFILE=""
OUTNAME=""

while getopts :t:p:r:1:2:o:q:h opt; do
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
		1)
			SEARCH_STRING=$OPTARG
			;;
		2)
			REPLACE_STRING=$OPTARG
			;;
		o)
			OUTNAME=$OPTARG
			;;
		q)
			FILE_PATTERN=$OPTARG
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
FASTQFILES=$(find $DATADIR -name "$FILE_PATTERN" -and -name "*${SEARCH_STRING}*") || exit

if [ "$FASTQFILES" == '' ]; then
	echo "Searching for files that match $FILE_PATTERN and ${SEARCH_STRING} in $DATADIR failed" >&2
	exit 1
fi

echo Specified files are: $FASTQFILES >&2
echo Specified mates are: ${FASTQFILES/$SEARCH_STRING/$REPLACE_STRING} >&2
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
	BASEFNAME=$(basename $F) || exit 1
	SECONDMATE=${F/$SEARCH_STRING/$REPLACE_STRING} || exit 1
	SECONDBASE=$(basename $SECONDMATE) || exit 1
	SAMOUT=${BASEFNAME%$SEARCH_STRING*}.sam || exit 1
	BAMOUT=${BASEFNAME%$SEARCH_STRING*}.bam || exit 1
	BAMS=$(echo $BAMS $BAMOUT) || exit 1
	if [ -e $BAMOUT ]; then samtools quickcheck $BAMOUT || rm $BAMOUT; fi #save time if execution was interrupted
	if [ ! -e $BAMOUT ]; then
		RGPU=$(head -n 1 $F | cut -d: -f3,4 --output-delimiter=.) || exit 1
		RGLB=$(expr $F : '.*\(M[0-9]*[abc]\)') || RGLB=$F || exit 1
		RGSM=$(expr $F : '.*\(M[0-9]*[abc]\)') || RGSM=$F || exit 1
		bwa mem -t ${CORES} -M -R '@RG\tID:'${RGSM}'\tPL:'${RGPL}'\tPU:'${RGPU}'\tLB:'${RGLB}'\tSM:'${RGSM} $REFERENCEFILE $F $SECONDMATE > $SAMOUT
		samtools sort -@ $CORES -o $BAMOUT -m 2G -T tmp $SAMOUT || exit 1
		rm $SAMOUT || exit 1
	fi
done

echo Merging . . . >&2
samtools merge -@ $CORES $OUTNAME $BAMS || exit 1

#Now clean up
rm $BAMS || exit

exit 0
