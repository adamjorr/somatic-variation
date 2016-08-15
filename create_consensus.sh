#!usr/bin/bash

USAGE="Usage: $0 [-t THREADS] [-f FILTER] [-c CHAINFILE] [-b BCFTOOLSFILE] [-o OUTPUT] reference.fasta in.bam"

CORES=48
FILTER='AC==1'
OUTPUT=consensus.fa
CHAINFILE=consensus.chain
BCFTOOLSFILE=bcftools_calls.vcf.gz

while getopts :t:h opt; do
	case $opt in
		t)
			CORES=$OPTARG
			;;
		b)
			BCFTOOLSFILE=$OPTARG
			;;
		c)
			CHAINFILE=$OPTARG
			;;
		f)
			FILTER=$OPTARG
			;;
		o)
			OUTPUT=$OPTARG
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

if [ $# -ne 2 ]; then			#if we forget arguments
	echo $USAGE	#remind us
	exit 1				#and exit with error
fi

REFFILE=$1
BAMFILE=$2

samtools mpileup -uf $REFFILE $BAMFILE | bcftools call --threads $CORES -mv -Ou | bcftools filter --threads $CORES -Oz -o $BCFTOOLSFILE -e $FILTER
tabix $BCFTOOLSFILE
cat $REFFILE | bcftools consensus $BCFTOOLSFILE -c $CHAINFILE > $OUTPUT


exit