#!/usr/bin/bash

#bwa_aligner.sh reference.fasta
#Takes reference and aligns all .fastq files in any subdirectory.
#This is a modified version of align_fastq which uses bwa instead of bowtie.
USAGE="Usage: $0 [-t THREADS] [-d TMPDIR] [-p RG_PLATFORM] [-q FILEPATTERN] [-1 FIRSTMATE] [-2 SECONDMATE] [-o outdir] -r reference.fa -i data/"



#Here are some things you might want to change:
RGPL=ILLUMINA #Default platform
FILE_PATTERN='*.fastq.gz' #pattern to find first set of reads
SEARCH_STRING='R1' #pattern for search/replace to find second set of reads
REPLACE_STRING='R2' #pattern for search/replace to substitute second set of reads
REFERENCEFILE=""
TMPOPTION=""
DATADIR=""
SENSITIVITY=""
OUTNAME=/dev/stdout

while getopts :t:d:p:r:1:2:o:i:q:s:h opt; do
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
			OUTDIR=$OPTARG
			;;
		i)
			DATADIR=$OPTARG
			;;
		q)
			FILE_PATTERN=$OPTARG
			;;
		s)
			SENSITIVITY=$OPTARG
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

if [ ! -d "$DATADIR" ]; then
	echo $USAGE >&2
	echo "$DATADIR is not a directory." >&2
	exit 1
fi


#Some variables
FASTQFILES=$(find $DATADIR -name "$FILE_PATTERN" -and -name "*${SEARCH_STRING}*")

if [ "$FASTQFILES" == '' ]; then
	echo "Searching for files that match $FILE_PATTERN and ${SEARCH_STRING} in $DATADIR failed" >&2
	exit 1
fi

if [ "$SENSITIVITY" != "" ]; then
	SENSITIVITY=$(echo -s $SENSITIVITY)
fi

echo Specified files are: $FASTQFILES >&2
echo Specified mates are: ${FASTQFILES/$SEARCH_STRING/$REPLACE_STRING} >&2
NUMFILES=( $FASTQFILES )
NUMFILES=${#NUMFILES[@]}

echo Making BAM files . . . >&2
#Make bamfiles from the FASTQs
PROGRESS=1
for F in $FASTQFILES; do
	echo Progress: $PROGRESS / $NUMFILES >&2
	BASEFNAME=$(basename $F)
	SECONDMATE=${F/$SEARCH_STRING/$REPLACE_STRING}
	SECONDBASE=$(basename $SECONDMATE)

	RGPU=$(zcat $F | head -n 1 | cut -d: -f3,4 --output-delimiter=.)
	RGLB=$(expr $F : '.*\(RL[0-9][0-9]\)') || RGLB=$F
	RGSM=$(expr $F : '.*\(RL[0-9][0-9]\)') || RGSM=$F

	qsub -v REF=$REFERENCEFILE,R1=$F,R2=$SECONDMATE,RGSM=$RGSM,RGLB=$RGLB,RGPU=$RGPU,BAMDIR=$OUTDIR HPC_jobs/ngm.job
	((PROGRESS++))
done

exit 0
