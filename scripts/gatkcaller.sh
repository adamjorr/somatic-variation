#!/usr/bin/bash

#http://tldp.org/LDP/abs/html/string-manipulation.html is a great guide for manipulating strings in bash
USAGE="Usage: $0 [-t THREADS] [-p PICARD_CMD] [-d TMPDIR] [-g  GATK_PATH] [-b BEDFILE] [-o OUTFILE] -r reference.fasta -i data.bam"

#entry point for this script is a deduplified bam.

#Here are some things you might want to change:
PICARD="picard" #How do I call picard on this system?
GATK=gatk #How do I call gatk on this system?
CORES=48
TMPOPTION=""
OUTFILE=/dev/stdout
NUMNS=30
BEDFILE=""
REFERENCEFILE=""
FILEIN=""
HET="0.025"
trap "exit 1" ERR

while getopts :t:p:g:d:b:o:r:i:h opt; do
	case $opt in
		t)
			CORES=$OPTARG
			;;
		p)
			PICARD=$OPTARG
			;;
		g)
			GATK=$OPTARG
			;;
		d)
			TMPOPTION=$OPTARG
			;;
		b)
			BEDFILE=$OPTARG
			;;
		o)
			OUTFILE=$OPTARG
			;;
		r)
			REFERENCEFILE=$OPTARG
			;;
		i)
			FILEIN=$OPTARG
			;;
		h)
			echo $USAGE >&2
			exit 1
			;;
	    \?)
			echo "Invalid option: -$OPTARG" >&2
			echo $USAGE >&2
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
	echo $USAGE >&2	#remind us
	exit 1				#and exit with error
fi

if [ "$REFERENCEFILE" == "" ]; then
	echo $USAGE >&2
	echo "Reference file required." >&2
	exit 1
fi

if [ "FILEIN" == "" ]; then
	echo $USAGE >&2
	echo "Input file required." >&2
	exit 1
fi

if [ "$TMPOPTION" != "" ]; then
	mkdir -p $TMPOPTION
fi

TMPDIR=$(mktemp -d --tmpdir=$TMPOPTION gatkcaller_tmp_XXXXXX)
trap "rm -rf $TMPDIR" EXIT INT TERM HUP

if [ "$BEDFILE" != "" ]; then
	BEDFILE=$(echo -XL $BEDFILE)
fi


##use GATK to do some analysis.
DEDUPLIFIEDBAM=${FILEIN}
REFERENCEDICT=${REFERENCEFILE%.*}.dict
SCATTEREDINTERVALDIR=$(mktemp --tmpdir=$TMPDIR scatteredIntervals_XXXXXX)
SCATTEREDFIRSTCALLDIR=$(mktemp --tmpdir=$TMPDIR scattered_first_calls_XXXXXX)
SUFFIXES=$(seq -f %02.0f 0 $((CORES-1)))
SCATTEREDFIRSTCALLS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i mktemp --tmpdir=$SCATTEREDFIRSTCALLDIR --suffix=.vcf.gz first_calls_{}_XXXXXX)
CMDFIRSTCALLS=$(echo $SCATTEREDFIRSTCALLS | tr ' ' '\n' | xargs -i echo I={})

JOINEDFIRSTCALLS=$(mktemp --tmpdir=$TMPDIR --suffix=.vcf joined_first_calls_XXXXXX)
SORTEDFIRSTCALLS=$(mktemp --tmpdir=$TMPDIR --suffix=.vcf sorted_first_calls_XXXXXX)
RECALIBRATEDBAM=$(mktemp --tmpdir=$TMPDIR --suffix=.bam recal_XXXXXX)
SCATTEREDOUTCALLDIR=$(mktemp --tmpdir=$TMPDIR scattered_output_calls_XXXXXX)
SCATTEREDOUTCALLS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i mktemp --tmpdir=$SCATTEREDOUTCALLDIR --suffix=.vcf out_call_{}_XXXXXX)
CMDOUTCALLS=$(echo $SCATTEREDOUTCALLS | tr ' ' '\n' | xargs -i echo I={})
RECALDATATABLE=$(mktemp --tmpdir=$TMPDIR --suffix=.table recal_data_XXXXXX)

if [ ! -e ${REFERENCEFILE}.fai ]; then
	samtools faidx $REFERENCEFILE
fi

if [ ! -e ${REFERENCEDICT} ]; then
	$PICARD CreateSequenceDictionary REFERENCE=${REFERENCEFILE} OUTPUT=${REFERENCEDICT}
fi

if [ ! -e ${DEDUPLIFIEDBAM%.bam}.bai ]; then
	$PICARD BuildBamIndex INPUT=${DEDUPLIFIEDBAM}
fi

#GATK TO RECALIBRATE QUAL SCORES + CALL VARIANTS
gatk SplitIntervals -R ${REFERENCEFILE} --scatter-count $CORES -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L scaffold_6 -L scaffold_7 -L scaffold_8 -L scaffold_9 -L scaffold_10 -L scaffold_11 -O ${SCATTEREDINTERVALDIR}
SCATTEREDINTERVALS=$(find ${SCATTEREDINTERVALDIR} -name '*.intervals')
parallel --halt 2 ${GATK} HaplotypeCaller --heterozygosity ${HET} -R ${REFERENCEFILE} -I $DEDUPLIFIEDBAM ${BEDFILE} -stand-call-conf 50 -ploidy 2 -L {1} -O {2} ::: ${SCATTEREDINTERVALS} :::+ ${SCATTEREDFIRSTCALLS}
$PICARD SortVcf ${CMDFIRSTCALLS} O=${SORTEDFIRSTCALLS} SEQUENCE_DICTIONARY=${REFERENCEDICT}
${GATK} BaseRecalibrator -I $DEDUPLIFIEDBAM -R ${REFERENCEFILE} ${BEDFILE} --known-sites $SORTEDFIRSTCALLS -O $RECALDATATABLE
${GATK} ApplyBQSR -I $DEDUPLIFIEDBAM -R ${REFERENCEFILE} --bqsr-recal-file $RECALDATATABLE -O $RECALIBRATEDBAM
parallel --halt 2 ${GATK} HaplotypeCaller -R ${REFERENCEFILE} -I $RECALIBRATEDBAM ${BEDFILE} --heterozygosity ${HET} -ploidy 2 -L {1} -O {2} ::: $SCATTEREDINTERVALS :::+ $SCATTEREDOUTCALLS
$PICARD SortVcf ${CMDOUTCALLS} O=${OUTFILE} SEQUENCE_DICTIONARY=${REFERENCEDICT}

exit 0
