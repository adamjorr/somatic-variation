#!usr/bin/bash

#stampy_realigner.sh reference.fasta data.bam out.bam
#Takes reference and realigns reads that are not in a proper pair
#or have quality lower than -q (default: 13) in the data.bam file.

USAGE="Usage: $0 [-t THREADS] [-d TMPDIR] [-q QUAL] [-o out.bam] -r reference.fasta -i data.bam"

CORES=48
TMPOPTION=""
QUAL=13
REFERENCEFILE=""
DATAFILE=""
OUTFILE=/dev/stdout

while getopts :t:d:q:r:i:o:h opt; do
	case $opt in
		t)
			CORES=$OPTARG
			;;
		d)
			TMPOPTION=$OPTARG
			;;
		q)
			QUAL=$OPTARG
			;;
		r)
			REFERENCEFILE=$OPTARG
			;;
		i)
			DATAFILE=$OPTARG
			;;
		o)
			OUTFILE=$OPTARG
			;;
		h)
			echo $USAGE >&2
			exit 1
			;;
	    \?)
			echo "Invalid option: -$OPTARG . Use -h for help." >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument. Use -h for help." >&2
			exit 1
			;;
	esac
done

shift $((OPTIND-1)) # get operands
if [ $# -ne 0 ]; then			#if we forget arguments
	echo $USAGE	>&2 #remind us
	exit 1				#and exit with error
fi

if [ "$REFERENCEFILE" == "" ]; then
	echo $USAGE >&2
	echo "Reference file required." >&2
	exit 1
fi

if [ "$DATAFILE" == "" ]; then
	echo $USAGE >&2
	echo "Input file required." >&2
	exit 1
fi

if [ ! -e "$DATAFILE" ]; then
	echo "$DATAFILE does not exist." >&2
	exit 1
fi

TMPDIR=$(mktemp -d --tmpdir=${TMPOPTION} stampy_realigner_tmp_XXXXXX)
SORTEDINFILE=$(mktemp --tmpdir=$TMPDIR --suffix=.bam sorted_in_XXX)
FIXEDINFILE=$(mktemp --tmpdir=$TMPDIR --suffix=.bam fixed_in_XXX)
MAPPEDREADS=$(mktemp --tmpdir=$TMPDIR --suffix=.bam mapped_XXX)
UNMAPPEDREADS=$(mktemp --tmpdir=$TMPDIR --suffix=.bam unmapped_XXX)

if [ ! -e ${REFERENCEFILE%.*}.stidx ]; then
	stampy -G ${REFERENCEFILE%.*} ${REFERENCEFILE} || exit 1
fi

if [ ! -e ${REFERENCEFILE%.*}.sthash ]; then
	stampy -g ${REFERENCEFILE%.*} -H ${REFERENCEFILE%.*} || exit 1
fi

echo SAMTOOLS PREPROCESSING >&2
samtools sort -@ ${CORES} -n -m 2G -T ${TMPDIR}/ -o $SORTEDINFILE $DATAFILE || exit 1
samtools fixmate $SORTEDINFILE $FIXEDINFILE || exit 1
rm $SORTEDINFILE || exit 1
samtools view -@ ${CORES} -b -h -q ${QUAL} -f 2 -o $MAPPEDREADS -U $UNMAPPEDREADS $FIXEDINFILE || exit 1
rm $FIXEDINFILE || exit 1



for GROUP in $(samtools view -H $UNMAPPEDREADS | grep ^@RG | cut -f2); do
	echo $GROUP >&2
	SANITARYGROUP=${GROUP//:/}
	SANITARYGROUP=${SANITARYGROUP//\//}
	SANITARYGROUP=${SANITARYGROUP//./}
	SANITARYGROUP=${SANITARYGROUP//\\/}
	GROUPFIFO=$(mktemp --tmpdir=$TMPDIR RG_${SANITARYGROUP}_XXX)
	mkfifo $GROUPFIFO
	FIFOS=$(echo $FIFOS $GROUPFIFO) || exit 1
    stampy -t ${CORES} -g ${REFERENCEFILE%.*} -h ${REFERENCEFILE%.*} -M $UNMAPPEDREADS --bamsortmemory=2000000000 --readgroup=${GROUP} |
    samtools view -@ ${CORES} -h -b -u | $GROUPFIFO
done

echo SAMTOOLS POSTPROCESSING >&2
MERGEDBAMS=$(mktemp --tmpdir=$TMPDIR --suffix=.bam merged_XXXXXX)
JOINEDMAPPEDREADS=$(mktemp --tmpdir=$TMPDIR --suffix=.bam joined_mapped_XXXXXX)


rm $UNMAPPEDREADS || exit 1
samtools merge -@ ${CORES} -n -c -p $MERGEDBAMS $FIFOS || exit 1
samtools merge -@ ${CORES} -n -c -p $OUTFILE $MERGEDBAMS $MAPPEDREADS || exit 1

rm -rf $TMPDIR

exit 0

