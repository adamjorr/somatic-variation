#!usr/bin/bash

#stampy_realigner.sh reference.fasta data.bam out.bam
#Takes reference and realigns all reads in the data.bam file.

USAGE="Usage: $0 [-t THREADS] reference.fasta data.bam out.bam"

CORES=48
REFERENCEFILE=$1
DATAFILE=$2
OUTFILE=$3

while getopts :t:h opt; do
	case $opt in
		t)
			CORES=$OPTARG
			;;
		h)
			echo $USAGE >&2
			exit 1
			;;
	    \?)
			echo "Invalid option: -$OPTARG ; use -h for help." >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument. Use -h for help." >&2
			exit 1
			;;
	esac
done

shift $((OPTIND-1)) # get operands
if [ $# -ne 3 ]; then			#if we forget arguments
	echo $USAGE	#remind us
	exit 1				#and exit with error
fi



if [ ! -e ${REFERENCEFILE%.*}.stidx ]; then
	stampy -G ${REFERENCEFILE%.*} ${REFERENCEFILE} || exit 1
fi

if [ ! -e ${REFERENCEFILE%.*}.sthash ]; then
	stampy -g ${REFERENCEFILE%.*} -H ${REFERENCEFILE%.*}
fi

for GROUP in $(samtools view -H $DATAFILE | grep ^@RG | cut -f2); do
	echo $GROUP
	SAMS=$(echo $SAMS ${GROUP#ID:}.sam) || exit
        stampy -t ${CORES} -g ${REFERENCEFILE%.*} -h ${REFERENCEFILE%.*} -M $DATAFILE --readgroup=${GROUP} > ${GROUP#ID:}.sam || exit 1
done

samtools merge -@ ${CORES} -n -c -p ${OUTFILE} $SAMS || exit 1
rm $SAMS
