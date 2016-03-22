#!usr/bin/bash

#bwa_aligner.sh reference.fasta
#Takes reference and aligns all .fastq files in any subdirectory.
#This is a modified version of align_fastq which uses bwa instead of bowtie.

if [ $# -ne 3 ]; then			#if we forget arguments
	echo "Usage: $0 reference.fasta data.bam out.bam"	#remind us
	exit 1				#and exit with error
fi

CORES=16
REFERENCEFILE=$1
DATAFILE=$2
OUTFILE=$3

for GROUP in $(samtools view -H fixmate_lowerthan13_unmapped.bam | grep ^@RG | cut -f2); do
		SAMS=$(echo $SAMS ${GROUP#ID:}.sam) || exit
        stampy -t 16 -g ../../grandis_reference/Egrandis_201 -h ../../grandis_reference/Egrandis_201 -M fixmate_lowerthan13_unmapped.bam --readgroup=${GROUP} > ${GROUP#ID:}.sam || exit 1
done

samtools merge -n -c -p $OUTFILE $SAMS
rm $SAMS