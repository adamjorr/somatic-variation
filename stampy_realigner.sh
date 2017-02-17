#!usr/bin/bash

#stampy_realigner.sh reference.fasta data.bam out.bam
#Takes reference and realigns reads that are not in a proper pair
#or have quality lower than -q (default: 13) in the data.bam file.

USAGE="Usage: $0 [-t THREADS] [-d TMPDIR] [-o out.bam] -r reference.fasta -i data.bam"

CORES=48
TMPOPTION=""
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
trap '[ -n "$(jobs -pr)" ] && kill $(jobs -pr)' INT TERM HUP
trap "exit 1" ERR
trap 'rm -rf $TMPDIR' EXIT


if [ ! -e ${REFERENCEFILE%.*}.stidx ]; then
	stampy -G ${REFERENCEFILE%.*} ${REFERENCEFILE}
fi

if [ ! -e ${REFERENCEFILE%.*}.sthash ]; then
	stampy -g ${REFERENCEFILE%.*} -H ${REFERENCEFILE%.*}
fi

READGROUPS=$(samtools view -H $DATAFILE | grep ^@RG | cut -f2)
# GROUPARRAY=( $READGROUPS )
# NUMOPERATIONS=$(( ${#GROUPARRAY[@]} + 1 ))
# CORESEACH=$(( $CORES / $NUMOPERATIONS ))
# if [ "$CORESEACH" -lt 1 ]; then
# 	CORESEACH=1
# fi

for GROUP in $READGROUPS; do
	echo $GROUP >&2
	SANITARYGROUP=${GROUP//:/}
	SANITARYGROUP=${SANITARYGROUP//\//}
	SANITARYGROUP=${SANITARYGROUP//./}
	SANITARYGROUP=${SANITARYGROUP//\\/}
	GROUPSAM=$(mktemp --suffix=.sam --tmpdir=$TMPDIR RG_${SANITARYGROUP}_XXX)
	SAMS=$(echo $SAMS $GROUPSAM)
    stampy -t ${CORES} -g ${REFERENCEFILE%.*} -h ${REFERENCEFILE%.*} --bamkeepgoodreads -M $DATAFILE --bamsortmemory=2000000000 --readgroup=${GROUP} |
    samtools sort -@ ${CORES} -T ${TMPDIR}/ -m 2G -O sam -o $GROUPSAM
done

samtools merge -@ ${CORES} -c -p $OUTFILE $FIFOS

exit 0

