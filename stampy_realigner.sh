#!usr/bin/bash

#stampy_realigner.sh reference.fasta data.bam out.bam
#Takes reference and realigns reads that are not in a proper pair
#or have quality lower than -q (default: 13) in the data.bam file.

USAGE="Usage: $0 [-t THREADS] [-d TMPDIR] [-q QUAL] reference.fasta data.bam out.bam"

CORES=48
TMPDIR=/tmp/
QUAL=13

while getopts :t:d:q:h opt; do
	case $opt in
		t)
			CORES=$OPTARG
			;;
		d)
			TMPDIR=$OPTARG
			;;
		q)
			QUAL=$OPTARG
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
REFERENCEFILE=$1
DATAFILE=$2
OUTFILE=$3



if [ ! -e ${REFERENCEFILE%.*}.stidx ]; then
	stampy -G ${REFERENCEFILE%.*} ${REFERENCEFILE} || exit 1
fi

if [ ! -e ${REFERENCEFILE%.*}.sthash ]; then
	stampy -g ${REFERENCEFILE%.*} -H ${REFERENCEFILE%.*} || exit 1
fi

echo SAMTOOLS PREPROCESSING >&2
samtools sort -@ ${CORES} -n -m 2G -o ${TMPDIR}/tmp.sorted.${DATAFILE} $DATAFILE || exit 1
samtools fixmate ${TMPDIR}/tmp.sorted.${DATAFILE} ${TMPDIR}/tmp.fixed.${DATAFILE} || exit 1
rm ${TMPDIR}/tmp.sorted.${DATAFILE} || exit 1
samtools view -@ ${CORES} -b -h -q ${QUAL} -f 2 -o ${TMPDIR}/tmp.mapped.bam -U ${TMPDIR}/tmp.unmapped.bam ${TMPDIR}/tmp.fixed.${DATAFILE} || exit 1
rm ${TMPDIR}/tmp.fixed.${DATAFILE} || exit 1

for GROUP in $(samtools view -H ${TMPDIR}/tmp.unmapped.bam | grep ^@RG | cut -f2); do
	echo $GROUP >&2
	BAMS=$(echo $BAMS ${TMPDIR}/${GROUP#ID:}.bam) || exit 1
        stampy -t ${CORES} -g ${REFERENCEFILE%.*} -h ${REFERENCEFILE%.*} -M ${TMPDIR}/tmp.unmapped.bam --readgroup=${GROUP} -o ${TMPDIR}/${GROUP#ID:}.stampy.sam || exit 1
        samtools sort -@ ${CORES} -n -m 2G -o ${TMPDIR}/${GROUP#ID:}.bam ${TMPDIR}/${GROUP#ID:}.stampy.sam || exit 1
        rm ${TMPDIR}/${GROUP#ID:}.stampy.sam || exit 1
done

echo SAMTOOLS POSTPROCESSING >&2

rm ${TMPDIR}/tmp.unmapped.bam || exit 1
samtools merge -@ ${CORES} -n -c -p ${TMPDIR}/tmp.stampy.${OUTFILE} $BAMS || exit 1
rm $BAMS || exit 1
samtools merge -@ ${CORES} -n -c -p ${TMPDIR}/tmp.namesorted.${OUTFILE} ${TMPDIR}/tmp.stampy.${OUTFILE} ${TMPDIR}/tmp.mapped.bam || exit 1
rm ${TMPDIR}/tmp.stampy.${OUTFILE} ${TMPDIR}/tmp.mapped.bam || exit 1
samtools sort -@ ${CORES} -m 2G -o ${OUTFILE} -O bam -T ${TMPDIR}/ ${TMPDIR}/tmp.namesorted.${OUTFILE} || exit 1
rm ${TMPDIR}/tmp.namesorted.${OUTFILE} || exit 1
