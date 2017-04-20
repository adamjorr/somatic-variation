#!/usr/bin/bash

USAGE="Usage: $0 [-t THREADS] [-d TMPDIR] [-f FILTER] [-c CHAINFILE] [-b BCFTOOLSFILE] -o OUTPUT -r reference.fasta -i in.bam"

CORES=48
FILTER='AC==1'
OUTPUT=""
BCFTOOLSFILE=""
BAMFILE=""
REFFILE=""
TMPOPTION=""

while getopts :t:d:b:c:f:o:i:r:h opt; do
	case $opt in
		t)
			CORES=$OPTARG
			;;
		d)
			TMPOPTION=$OPTARG
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
		i)
			BAMFILE=$OPTARG
			;;
		r)
			REFFILE=$OPTARG
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

if [ $# -ne 0 ]; then			#if we forget arguments
	echo $USAGE	>&2 #remind us
	exit 1				#and exit with error
fi

if [ "$BAMFILE" == "" ]; then
	echo $USAGE >&2
	echo "Input file required."
	exit 1
fi

if [ "$OUTPUT" == "" ]; then
	echo $USAGE >&2
	echo "Output file required."
	exit 1
fi

if [ "$CHAINFILE" == "" ]; then
	CHAINFILE=$(basename ${OUTPUT}.chain)
fi

if [ "$REFFILE" == "" ]; then
	echo $USAGE >&2
	echo "Reference file required."
	exit 1
fi

TMPDIR=$(mktemp -d --tmpdir=$TMPOPTION $(basename $0)_tmp_XXXXXX)
trap "rm -rf $TMPDIR" EXIT INT TERM HUP
trap "exit 1" ERR

if [ "$BCFTOOLSFILE" == "" ]; then
	BCFTOOLSFILE=$(mktemp --tmpdir=$TMPDIR --suffix=.vcf.gz bcftools_calls_XXXXXX)
fi

CHRS=$(samtools view -H $BAMFILE | grep ^@SQ | cut -f 2 | sed -n -e 's/^.*://p')
CHRSARY=( $CHRS )
NUMCHR=${#CHRSARY[@]}
if [ $CORES -gt $NUMCHR ]; then
	CORES=$NUMCHR
fi

CHRFILES=$(echo ${CHRS[@]} | xargs -n 1 | xargs -n 1 -I{} mktemp --tmpdir=$TMPDIR --suffix=.vcf.gz {}_XXXXXX)
DEDUPLIFIEDBAM=$(mktemp --tmpdir=$TMPDIR --suffix=.bam dedup_XXX)
METRICFILE=$(mktemp --tmpdir=$TMPDIR --suffix=.txt metrics_XXX)
$PICARD MarkDuplicates INPUT=${BAMFILE} OUTPUT=${DEDUPLIFIEDBAM} METRICS_FILE=$METRICFILE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
samtools index ${DEDUPLIFIEDBAM}

export DEDUPLIFIEDBAM
export REFFILE
export FILTER
parallel --link -j $CORES --halt now,fail=1 --env DEDUPLIFIEDBAM --env REFFILE --env FILTER \
'samtools mpileup -r {1} -guf ${REFFILE} ${DEDUPLIFIEDBAM} | bcftools call -mv -Ou | bcftools filter -Oz -o {2} -e ${FILTER}' \
::: ${CHRS[@]} ::: ${CHRFILES}

bcftools concat -Oz -o $BCFTOOLSFILE $CHRFILES
tabix $BCFTOOLSFILE
cat $REFFILE | bcftools consensus $BCFTOOLSFILE -c $CHAINFILE > $OUTPUT

exit 0
