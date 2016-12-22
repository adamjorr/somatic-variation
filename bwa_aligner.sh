#!usr/bin/bash

#bwa_aligner.sh reference.fasta
#Takes reference and aligns all .fastq files in any subdirectory.
#This is a modified version of align_fastq which uses bwa instead of bowtie.
USAGE="Usage: $0 [-t THREADS] [-d TMPDIR] [-p RG_PLATFORM] [-q FILEPATTERN] [-1 FIRSTMATE] [-2 SECONDMATE] [-o out.bam] -r reference.fa -i data/"




#Here are some things you might want to change:
RGPL=ILLUMINA #We assume Illumina; if we're wrong, change it here.
CORES=48
FILE_PATTERN='*.fastq' #pattern to find first set of reads
SEARCH_STRING='R1' #pattern for search/replace to find second set of reads
REPLACE_STRING='R2' #pattern for search/replace to substitute second set of reads
REFERENCEFILE=""
TMPOPTION=""
DATADIR=""
OUTNAME=/dev/stdout

while getopts :t:p:r:1:2:o:i:q:h opt; do
	case $opt in
		t)
			CORES=$OPTARG
			;;
		d)
			TMPOPTION=$OPTARG
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
		i)
			DATADIR=$OPTARG
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
if [ $# -ne 0 ]; then			#if we forget arguments
	echo $USAGE	>&2#remind us
	exit 1				#and exit with error
fi

if [ "$REFERENCEFILE" == "" ]; then
	echo $USAGE >&2
	echo "Reference file required." >&2
	exit 1
fi

if [ "$DATADIR" == "" ]; then
	echo $USAGE >&2
	echo "Input directory required." >&2
	exit 1
fi

if [ ! -d "$DATADIR"]; then
	echo $USAGE >&2
	echo "$DATADIR is not a directory." >&2
	exit 1
fi

TMPDIR=$(mktemp -d --tmpdir=$TMPOPTION $0_tmp_XXXXXX)
trap "rm -rf $TMPDIR" EXIT INT TERM HUP
trap "exit 1" ERR


#Some variables
FASTQFILES=$(find $DATADIR -name "$FILE_PATTERN" -and -name "*${SEARCH_STRING}*")

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
	BASEFNAME=$(basename $F)
	SECONDMATE=${F/$SEARCH_STRING/$REPLACE_STRING}
	SECONDBASE=$(basename $SECONDMATE)
	BAMOUT=$(mktemp --tmpdir=$TMPDIR --suffix=.bam ${BASEFNAME%$SEARCH_STRING*}_XXXXXXXXXX)
	BAMS=$(echo $BAMS $BAMOUT)
	RGPU=$(head -n 1 $F | cut -d: -f3,4 --output-delimiter=.)
	RGLB=$(expr $F : '.*\(M[0-9]*[abc]\)') || RGLB=$F
	RGSM=$(expr $F : '.*\(M[0-9]*[abc]\)') || RGSM=$F
	bwa mem -t ${CORES} -M -R '@RG\tID:'${RGSM}'\tPL:'${RGPL}'\tPU:'${RGPU}'\tLB:'${RGLB}'\tSM:'${RGSM} $REFERENCEFILE $F $SECONDMATE | 
		samtools sort -@ $CORES -o $BAMOUT -O bam -m 2G -T ${TMPDIR}/
done

echo Merging . . . >&2
samtools merge -@ $CORES $OUTNAME $BAMS

exit 0
