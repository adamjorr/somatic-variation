#!/usr/bin/bash
# NCI version of gatkcaller.sh
# It submits it as a job to the cluster

#http://tldp.org/LDP/abs/html/string-manipulation.html is a great guide for manipulating strings in bash
USAGE="Usage: $0 [-t THREADS] [-p PICARD_CMD] [-d TMPDIR] [-g  GATK_PATH] [-b BEDFILE] [-o OUTFILE] -r reference.fasta -i deduplified.bam"

#Here are some things you might want to change:
module load samtools/1.4
module load java
JOBS=/short/xf1/src_big/somatic-variation/HPC_jobs
PICARD="java -jar $HOME/apps/bin/picard.jar" #How do I call picard on this system?
GATK=/short/xf1/src_big/GATK/GenomeAnalysisTK.jar #Location of your GATK jar
CORES=48		#NCI submitted jobs will use 16 cores in the job script
TMPOPTION=""
OUTFILE=/dev/stdout
NUMNS=30
BEDFILE=""
REFERENCEFILE=""
FILEIN=""
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
			DEDUPLIFIEDBAM=$OPTARG
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


TMPDIR=$(mktemp -d --tmpdir=$TMPOPTION gatkcaller_tmp_XXXXXX)
# trap "rm -rf $TMPDIR" EXIT INT TERM HUP

if [ "$BEDFILE" != "" ]; then
	BEDFILE=$(echo -XL $BEDFILE)
fi


### Set up a crap-ton of temporary directories and files!
REFERENCEDICT=${REFERENCEFILE%.*}.dict
FULLINTERVALS=$(mktemp --tmpdir=$TMPDIR --suffix=.interval_list fullIntervals_XXX)
SCATTEREDINTERVALDIR=$(mktemp -d --tmpdir=$TMPDIR scatteredIntervals_XXXXXX)
SCATTEREDFIRSTCALLDIR=$(mktemp -d --tmpdir=$TMPDIR scattered_first_calls_XXX)
SUFFIXES=$(seq -f %02.0f 0 $((CORES-1)))
SCATTEREDFIRSTCALLS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i mktemp --tmpdir=$SCATTEREDFIRSTCALLDIR --suffix=.vcf.gz first_calls_{}_XXXXXX)

CMDFIRSTCALLS=$(echo $SCATTEREDFIRSTCALLS | tr ' ' '\n' | xargs -i echo I={})
JOINEDFIRSTCALLS=$(mktemp --tmpdir=$TMPDIR --suffix=.vcf joined_first_calls_XXX)
SORTEDFIRSTCALLS=$(mktemp --tmpdir=$TMPDIR --suffix=.vcf sorted_first_calls_XXX)
RECALIBRATEDBAM=$(mktemp --tmpdir=$TMPDIR --suffix=.bam recal_XXX)
SCATTEREDOUTCALLDIR=$(mktemp -d --tmpdir=$TMPDIR scattered_output_calls_XXX)
SCATTEREDOUTCALLS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i mktemp --tmpdir=$SCATTEREDOUTCALLDIR --suffix=.vcf.gz out_call_{}_XXXXXX)
CMDOUTCALLS=$(echo $SCATTEREDOUTCALLS | tr ' ' '\n' | xargs -i echo I={})
OUTCALLS=$(mktemp --tmpdir=$TMPDIR --suffix=.vcf out_calls_XXX)
RECALDATATABLE=$(mktemp --tmpdir=$TMPDIR --suffix=.table recal_data_XXX)


if [ ! -e ${REFERENCEFILE}.fai ]; then
	samtools faidx $REFERENCEFILE			# quick
fi

if [ ! -e ${REFERENCEDICT} ]; then
	echo "$PICARD CreateSequenceDictionary REFERENCE=${REFERENCEFILE} OUTPUT=${REFERENCEDICT}"
	$PICARD CreateSequenceDictionary REFERENCE=${REFERENCEFILE} OUTPUT=${REFERENCEDICT}	# very quick
fi

#GATK TO RECALIBRATE QUAL SCORES + CALL VARIANTS
echo "$PICARD ScatterIntervalsByNs R=${REFERENCEFILE} OT=ACGT MAX_TO_MERGE=${NUMNS} O=${FULLINTERVALS}"
$PICARD ScatterIntervalsByNs R=${REFERENCEFILE} OT=ACGT MAX_TO_MERGE=${NUMNS} O=${FULLINTERVALS}	#very quick
$PICARD IntervalListTools I=${FULLINTERVALS} SCATTER_COUNT=$CORES O=${SCATTEREDINTERVALDIR}		#very quick. For each core it makes a tmp directory containing a subset of intervals.
													# HOWEVER, this is a bottleneck on the head node. Can't execute with CORES > 48 or so. We need to job it :( 

SCATTEREDINTERVALS=$(find ${SCATTEREDINTERVALDIR} -name '*.interval_list')

#  run GATK on a whole lot of intervals using parallel. Each call to GATK uses one ScatteredInterval as input and its corresponding ScatterFirstcall as output
# DK: this *could* be massively multiplexed on NCI. We will start by just submitting the parallel call to a node with lots of CPUs and then see if it is worth optimizing further.

echo "Time to qsub GATK.job"

ARRAY_SI=($SCATTEREDINTERVALS)
ARRAY_SFC=($SCATTEREDFIRSTCALLS)

for ((i=0; i<$CORES; i++));
do
	SI=${ARRAY_SI[$i]}
	SFC=${ARRAY_SFC[$i]}
	echo "qsub GATK.job "${i}
	qsub -v REF=$REFERENCEFILE,BAM=$DEDUPLIFIEDBAM,BED=$BEDFILE,SCATTEREDINTERVAL=$SI,SCATTEREDCALL=$SFC $JOBS/GATK.job
	sleep 0.1
done


#parallel --halt 2 java -jar ${GATK} -T HaplotypeCaller --pair_hmm_implementation LOGLESS_CACHING -R ${REFERENCEFILE} -I $DEDUPLIFIEDBAM -L {1} ${BEDFILE} \
# -stand_call_conf 50 -ploidy 2 -o {2} ::: $SCATTEREDINTERVALS :::+ $SCATTEREDFIRSTCALLS

# java -cp ${GATK} org.broadinstitute.gatk.tools.CatVariants -R ${REFERENCEFILE} --outputFile ${JOINEDFIRSTCALLS} ${CMDFIRSTCALLS}
# bcftools concat -a -Ov -o ${JOINEDFIRSTCALLS} ${SCATTEREDFIRSTCALLS}

: << 'BLOCKCOMMENT'

$PICARD SortVcf ${CMDFIRSTCALLS} O=${SORTEDFIRSTCALLS} SEQUENCE_DICTIONARY=${REFERENCEDICT}
rm $SCATTEREDFIRSTCALLS ${SORTEDFIRSTCALLS}.idx
# rm $JOINEDFIRSTCALLS

java -jar ${GATK} -T BaseRecalibrator -nct $CORES -I $DEDUPLIFIEDBAM -R ${REFERENCEFILE} ${BEDFILE} --knownSites $SORTEDFIRSTCALLS -o $RECALDATATABLE
rm $SORTEDFIRSTCALLS

java -jar ${GATK} -T PrintReads -nct $CORES -I $DEDUPLIFIEDBAM -R ${REFERENCEFILE} -BQSR $RECALDATATABLE -EOQ -o $RECALIBRATEDBAM
rm $DEDUPLIFIEDBAM $RECALDATATABLE

ARRAY_SOC=($SCATTEREDOUTCALLS)

for ((i=0; i<$CORES; i++));
do
        SI=${ARRAY_SI[$i]}
        SFC=${ARRAY_SFC[$i]}
        echo "qsub GATK.job "${i}
        qsub -v REF=$REFERENCEFILE,BAM=$DEDUPLIFIEDBAM,BED=$BEDFILE,SCATTEREDINTERVAL=$SI,SCATTEREDCALL=$SOC $JOBS/GATK.job
        sleep 0.1
done


#parallel --halt 2 java -jar ${GATK} -T HaplotypeCaller --pair_hmm_implementation LOGLESS_CACHING -R ${REFERENCEFILE} -I $RECALIBRATEDBAM -L {1} ${BEDFILE} -ploidy 2 -o {2} ::: $SCATTEREDINTERVALS :::+ $SCATTEREDOUTCALLS

rm $SCATTEREDINTERVALS $RECALIBRATEDBAM
# java -cp ${GATK} org.broadinstitute.gatk.tools.CatVariants -R ${REFERENCEFILE} --outputFile ${OUTFILE} ${CMDOUTCALLS}
# bcftools concat -a -Ov -o ${OUTCALLS} ${SCATTEREDOUTCALLS}

$PICARD SortVcf ${CMDOUTCALLS} O=${OUTFILE} SEQUENCE_DICTIONARY=${REFERENCEDICT}

BLOCKCOMMENT

exit 0
