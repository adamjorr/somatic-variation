#!usr/bin/bash

#stampy_realigner.sh reference.fasta data.bam out.bam
#Takes reference and realigns all reads in the data.bam file.

if [ $# -ne 3 ]; then			#if we forget arguments
	echo "Usage: $0 reference.fasta data.bam out.bam"	#remind us
	exit 1				#and exit with error
fi

CORES=16
REFERENCEFILE=$1
DATAFILE=$2
OUTFILE=$3

if ! [ -e ${REFERENCEFILE%%.*}.stidx ]
	stampy -G ${REFERENCEFILE%%.*} ${REFERENCEFILE} || exit 1
fi

if ! [ -e ${REFERENCEFILE%%.*}.sthash ]
	stampy -g ${REFERENCEFILE%%.*} -H ${REFERENCEFILE%%.*}
fi

for GROUP in $(samtools view -H fixmate_lowerthan13_unmapped.bam | grep ^@RG | cut -f2); do
		SAMS=$(echo $SAMS ${GROUP#ID:}.sam) || exit
        stampy -t 16 -g ${REFERENCEFILE%%.*} -h ${REFERENCEFILE%%.*} -M $DATAFILE --readgroup=${GROUP} > ${GROUP#ID:}.sam || exit 1
done

samtools merge -@ 16 -n -c -p --O BAM ${OUTFILE} $SAMS
rm $SAMS