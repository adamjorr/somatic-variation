#!/usr/bin/env bash
#bash clean_reads.sh -t THREADS -d DEST_DIRECTORY -k KMER_SIZE -f FILE_PATTERN -s SEARCH_STRING -r REPLACE_STRING [-m MAX_MEMORY] [-c COVERAGE] READ_DIRECTORY 
USAGE="Usage: $0 [-t THREADS] [-d DEST_DIRECTORY] [-k KMER_SIZE] [-f FILE_PATTERN] [-s SEARCH_STRING] [-r REPLACE_STRING] [-m MAX_MEMORY] [-c COVERAGE] READ_DIRECTORY"

RCORRECTOR="perl run_rcorrector.pl"
LOAD_COUNTING="scripts/load-into-counting.py"
SLICE_BY_COV="scripts/slice-paired-reads-by-coverage.py"
THREADS=48
DEST_DIRECTORY=./cleaned_reads
KMER_SIZE=32
FILE_PATTERN='*R1*.fastq' #pattern to find first set of reads
SEARCH_STRING='R1' #pattern for search/replace to find second set of reads
REPLACE_STRING='R2' #pattern for search/replace to find second set of reads
MAX_MEMORY=64e9 #max memory to be given to khmer
COVERAGE=40000 #max coverage tolerable  (see khmer slice-reads-by-coverage)

while getopts :t:d:k:f:s:r:m:c:h opt; do
	case $opt in
		t)
			echo "t was set to $OPTARG" >&2
			THREADS=$OPTARG
			;;
		d)
			echo "d was set to $OPTARG" >&2
			DEST_DIRECTORY=$OPTARG
			;;
		k)
			echo "k was set to $OPTARG" >&2
			KMER_SIZE=$OPTARG
			;;
		f)
			echo "f was set to $OPTARG" >&2
			FILE_PATTERN=$OPTARG
			;;
		s)
			echo "s was set to $OPTARG" >&2
			SEARCH_STRING=$OPTARG
			;;
		r)
			echo "r was set to $OPTARG" >&2
			REPLACE_STRING=$OPTARG
			;;
		m)
			echo "m was set to $OPTARG" >&2
			MAX_MEMORY=$OPTARG
			;;
		c)
			echo "c was set to $OPTARG" >&2
			COVERAGE=$OPTARG
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

if [ $# -ne 1 ]; then
	echo $USAGE >&2
	exit 1
fi

if [ ! -d $1 ] ; then
	echo "Invalid argument '$1': Not a directory or unreadable."
	exit 1
fi

READ_DIRECTORY=$1
mkdir -p $DEST_DIRECTORY
mkdir -p ${DEST_DIRECTORY}/corrected
mkdir -p ${DEST_DIRECTORY}/sliced

# Run Rcorrector with the script included in the program:
for F in $(find $DEST_DIRECTORY -name $FILE_PATTERN); do
	$RCORRECTOR -1 $F -2 ${F/$SEARCH_STRING/$REPLACE_STRING} -k $KMER_SIZE -t $THREADS -od ${DEST_DIRECTORY}/corrected/ || exit 1
done

# Build a graph and filter on estimated coverage using Khmer. Check [this fork](https://github.com/adamjorr/khmer):
$LOAD_COUNTING -ksize $KMER_SIZE -T $THREADS -M 64e9 khmer_count.graph ${DEST_DIRECTORY}/corrected/*.fq || exit 1

# For a parallelized version, use:
parallel -j $THREADS 'F={}; G={/.}; $SLICE_BY_COV -M $COVERAGE khmer_count.graph $F ${F/$SEARCH_STRING/$REPLACE_STRING} ${DEST_DIRECTORY}/sliced/${G}_sliced.fq ${DEST_DIRECTORY}/sliced/${G/$SEARCH_STRING/$REPLACE_STRING}_sliced.fq ${DEST_DIRECTORY}/sliced/${G/$SEARCH_STRING/}_singletons.fq' ::: $(find $DEST_DIRECTORY/corrected/ -name "*${SEARCH_STRING}*.fq") || exit 1

exit
