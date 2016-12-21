#!usr/bin/bash
#bash vcf2fa.sh <file.vcf >out.fa
#Takes vcf from STDIN (or -i) and prints a fasta file to STDOUT (or -o)

USAGE="$0 [-d tmp_directory/] [-i file.vcf] [-o out.fa]"

INFILE="-"
OUTFILE="/dev/stdout"
trap "exit 1" ERR
trap '[ -n "$(jobs -pr)" ] && kill $(jobs -pr); rm -rf $TMPDIR' EXIT INT TERM HUP

while getopts d:i:o:h opt; do
	case $opt in
		d)
			TMPOPT=$OPTARG
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

$TMPDIR=$(mktemp --tmpdir=$TMPOPT vcf2fa_tmp_XXXXXXXX)

$FIFONAME=$(mktemp -u --tmpdir=$TMPDIR --suffix=.tab $0_pipe_XXXXXXXX)
mkfifo $FIFONAME;

cat $INFILE | vcf-to-tab >$FIFONAME &
perl vcf_tab_to_fasta_alignment.pl -i $FIFONAME >$OUTFILE
