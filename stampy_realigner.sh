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

trap '[ -n "$(jobs -pr)" ] && kill $(jobs -pr); rm -rf $TMPDIR' EXIT INT TERM HUP
trap "exit 1" ERR


if [ ! -e ${REFERENCEFILE%.*}.stidx ]; then
	stampy -G ${REFERENCEFILE%.*} ${REFERENCEFILE}
fi

if [ ! -e ${REFERENCEFILE%.*}.sthash ]; then
	stampy -g ${REFERENCEFILE%.*} -H ${REFERENCEFILE%.*}
fi

echo SAMTOOLS PREPROCESSING >&2
samtools sort -@ ${CORES} -n -m 2G -T ${TMPDIR}/ -o $SORTEDINFILE $DATAFILE
samtools fixmate $SORTEDINFILE $FIXEDINFILE
rm $SORTEDINFILE
samtools view -@ ${CORES} -b -h -q ${QUAL} -f 2 -f 4 -o $MAPPEDREADS -U $UNMAPPEDREADS $FIXEDINFILE
rm $FIXEDINFILE



for GROUP in $(samtools view -H $UNMAPPEDREADS | grep ^@RG | cut -f2); do
	echo $GROUP >&2
	SANITARYGROUP=${GROUP//:/}
	SANITARYGROUP=${SANITARYGROUP//\//}
	SANITARYGROUP=${SANITARYGROUP//./}
	SANITARYGROUP=${SANITARYGROUP//\\/}
	GROUPFIFO=$(mktemp -u --suffix=.bam --tmpdir=$TMPDIR RG_${SANITARYGROUP}_XXX)
	mkfifo $GROUPFIFO
	FIFOS=$(echo $FIFOS $GROUPFIFO)
    stampy -t ${CORES} -g ${REFERENCEFILE%.*} -h ${REFERENCEFILE%.*} -M $UNMAPPEDREADS --bamsortmemory=2000000000 --readgroup=${GROUP} |
    samtools sort -@ ${CORES} -T ${TMPDIR}/ -m 2G -n -O bam > $GROUPFIFO &
done

samtools merge -@ ${CORES} -n -c -p $OUTFILE $MAPPEDREADS $FIFOS

exit 0

