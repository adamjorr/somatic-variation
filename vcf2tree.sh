#!usr/bin/bash
#bash vcf2tree.sh file.vcf out.pdf
#Takes file.vcf, filters it using replicate info at various stringencies, and plots the trees.

USAGE="$0 [-t THREADS] -g GROUPBY -i file.vcf -o tree.pdf"

THREADS=4
INFILE=""
OUTFILE=""
GROUPBY=""
trap "exit 1" ERR

while getopts t:g:i:o:h opt; do
	case $opt in
		t)
			THREADS=$OPTARG
			;;
		g)
			GROUPBY=$OPTARG
			;;
		i)
			INFILE=$OPTARG
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

if [ "$INFILE" == "" ]; then
	echo $USAGE >&2
	echo "Input file required." >&2
	exit 1
fi

if [ "$GROUPBY" == "" ]; then
	echo $USAGE >&2
	echo "Size of sample groups required." >&2
	exit 1
fi

if [ "$OUTFILE" == "" ]; then
	echo $USAGE >&2
	echo "Out file required." >&2
	exit 1
fi

mkdir -p vcf_processing
cd vcf_processing

perl filt_with_replicates.pl -s -g $GROUPBY <$INFILE >filtered.vcf
vcf-to-tab <filtered.vcf >filtered.tab
perl vcf_tab_to_fasta_alignment.pl -i filtered.tab > cleaned.fasta
rm filtered.tab_clean
python2 diploidify.py -i cleaned.fasta -t fasta -o cleaned.dip.phylip-relaxed -p phylip-relaxed -v

mkdir -p tree
cd tree

raxmlHPC -T $THREADS -f a -s ../cleaned.dip.phylip-relaxed -n nwk -m ASC_GTRGAMMA --asc-corr=lewis -p 12345 -x 12345 -# 100

cd ..

Rscript plot_tree.R tree/RAxML_bestTree.nwk $OUTFILE

cd ..


exit
