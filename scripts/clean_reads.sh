#!/usr/bin/env bash
#bash clean_reads.sh -t THREADS -d DEST_DIRECTORY -k KMER_SIZE -f FILE_SUFFIX -s SEARCH_STRING -r REPLACE_STRING [-m MAX_MEMORY] [-c COVERAGE] READ_DIRECTORY 
USAGE="Usage: $0 [-t THREADS] [-d DEST_DIRECTORY] [-k KMER_SIZE] [-f FILE_SUFFIX] [-1 SEARCH_STRING] [-2 REPLACE_STRING] [-m MAX_MEMORY] [-c COVERAGE] -i READ_DIRECTORY/"

RCORRECTOR="run_rcorrector.pl"
LOAD_COUNTING="load-into-counting.py"
SLICE_BY_COV="slice-paired-reads-by-coverage.py"
THREADS=48
DEST_DIRECTORY=./cleaned_reads
KMER_SIZE=32
FILE_SUFFIX='.fastq' #pattern to find first set of reads
SEARCH_STRING='R1' #pattern for search/replace to find second set of reads
REPLACE_STRING='R2' #pattern for search/replace to substitute second set of reads
MAX_MEMORY=64e9 #max memory to be given to khmer
COVERAGE=40000 #max coverage tolerable  (see khmer slice-reads-by-coverage)
READ_DIRECTORY=""
SLICE_THREADS=""

while getopts :t:d:k:f:1:2:m:c:i:s:h opt; do
	case $opt in
		t)
			THREADS=$OPTARG
			;;
		d)
			DEST_DIRECTORY=$OPTARG
			;;
		k)
			KMER_SIZE=$OPTARG
			;;
		f)
			FILE_SUFFIX=$OPTARG
			;;
		1)
			SEARCH_STRING=$OPTARG
			;;
		2)
			REPLACE_STRING=$OPTARG
			;;
		m)
			MAX_MEMORY=$OPTARG
			;;
		c)
			COVERAGE=$OPTARG
			;;
		i)
			READ_DIRECTORY=$OPTARG
			;;
		s)
			SLICE_THREADS=$OPTARG
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

shift $((OPTIND-1))

if [ $# -ne 0 ]; then
	echo $USAGE >&2
	exit 1
fi

if [ "$READ_DIRECTORY" == "" ]; then
	echo $USAGE >&2
	echo "Input directory required." >&2
	exit 1
fi

if [ ! -d $READ_DIRECTORY ] ; then
	echo $USAGE >&2
	echo "Invalid argument ${READ_DIRECTORY}: Not a directory or unreadable." >&2
	exit 1
fi

if [ "$SLICE_THREADS" == "" ]; then
	SLICE_THREADS=$((THREADS/6))
fi

if [ $SLICE_THREADS -eq 0 ]; then
	SLICE_THREADS=1
fi

trap "exit 1" ERR

CORR_SUFFIX=.cor${FILE_SUFFIX//.fastq/.fq}
mkdir -p $DEST_DIRECTORY
mkdir -p ${DEST_DIRECTORY}/corrected
mkdir -p ${DEST_DIRECTORY}/sliced

# Run Rcorrector with the script included in the program:
for F in $(find $READ_DIRECTORY -name "*$FILE_SUFFIX" -and -name "*$SEARCH_STRING*"); do
	$RCORRECTOR -1 $F -2 ${F/$SEARCH_STRING/$REPLACE_STRING} -k $KMER_SIZE -t $THREADS -od ${DEST_DIRECTORY}/corrected/
done

# Build a graph with khmer
$LOAD_COUNTING -k $KMER_SIZE -T $THREADS -M $MAX_MEMORY ${DEST_DIRECTORY}/khmer_count.graph ${DEST_DIRECTORY}/corrected/*${CORR_SUFFIX}

# filter on $COVERAGE in parallel using $SLICE_BY_COV
export SLICE_BY_COV
export COVERAGE
export SEARCH_STRING
export REPLACE_STRING
export DEST_DIRECTORY
export FILE_SUFFIX
export CORR_SUFFIX
parallel -j $SLICE_THREADS --env SLICE_BY_COV --env COVERAGE --env SEARCH_STRING --env REPLACE_STRING \
--env DEST_DIRECTORY --env FILE_SUFFIX --env CORR_SUFFIX 'F={}; G={/.}; \
$SLICE_BY_COV -M $COVERAGE ${DEST_DIRECTORY}/khmer_count.graph $F ${F/$SEARCH_STRING/$REPLACE_STRING} \
${DEST_DIRECTORY}/sliced/${G}_sliced${FILE_SUFFIX} \
${DEST_DIRECTORY}/sliced/${G/$SEARCH_STRING/$REPLACE_STRING}_sliced${FILE_SUFFIX} \
${DEST_DIRECTORY}/sliced/${G/$SEARCH_STRING/}_singletons${FILE_SUFFIX}' \
::: $(find $DEST_DIRECTORY/corrected/ -name "*$CORR_SUFFIX" -and -name "*${SEARCH_STRING}*")

exit
