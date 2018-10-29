#!/usr/bin/bash
#bash vcf2tree.sh file.vcf out.pdf
#Takes a VCF, filters it using replicate info, and constructs a tree with RAxML

USAGE="$0 [-t THREADS] [-r raxmlHPC] [-i file.vcf] [-o tree.nwk] -g GROUPBY"

THREADS=4
INFILE="-"
OUTFILE="/dev/stdout"
GROUPBY=""
RAXMLCALL=raxmlHPC
trap "exit 1" ERR
trap 'rm -rf $TMPDIR' EXIT INT TERM HUP

while getopts t:g:r:i:o:h opt; do
	case $opt in
		t)
			THREADS=$OPTARG
			;;
		g)
			GROUPBY=$OPTARG
			;;
		r)
			RAXMLCALL=$OPTARG
			;;
		i)
			INFILE=$OPTARG
			;;
		d)
			TMPOPT=$OPTARG
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

if [ "$GROUPBY" == "" ]; then
	echo $USAGE >&2
	echo "Size of sample groups required." >&2
	exit 1
fi

TMPDIR=$(mktemp -d --tmpdir=$TMPOPT vcf2tree_tmp_XXXXXXXX)
DIPLOIDIFIED=$(mktemp --tmpdir=$TMPDIR --suffix=.fa diploidified_tmp_XXXXXX)
TREEFOLDER=$(mktemp -d --tmpdir=$TMPDIR vcf2tree_tree_tmp_XXXXXX)

cat $INFILE | filt_with_replicates.pl -s -g $GROUPBY | vcf2fa.sh | diploidify.py -v > $DIPLOIDIFIED
${RAXMLCALL} -T $THREADS -f a -s $DIPLOIDIFIED -n nwk -m ASC_GTRGAMMA -w $TREEFOLDER --asc-corr=lewis -p 12345 -x 12345 -# 100
cat ${TREEFOLDER}/RAxML_bestTree.nwk >$OUTFILE

exit 0
