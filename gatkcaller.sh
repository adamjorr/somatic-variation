#!usr/bin/bash
##TODO: REMOVE ALIGNMENT COMPONENTS, UPDATE DOCS + VARIABLES
#aligntoref.sh reference.fasta data.bam
#Takes reference and aligns all .fastq files in any subdirectory and calls SNPs with GATK.
#requires samtools, picard-tools, the GATK, and bowtie2

#http://tldp.org/LDP/abs/html/string-manipulation.html is a great guide for manipulating strings in bash
USAGE="Usage: $0 [-t THREADS] [-p PICARD_CMD] [-d TMPDIR] [-g  GATK_PATH] [-o OUTFILE] reference.fasta data.bam"

#Here are some things you might want to change:
PICARD="picard" #How do I call picard on this system?
GATK=~/bin/GenomeAnalysisTK.jar #Location of your GATK jar
CORES=48
TMPOPTION=""
OUTFILE=/dev/stdout
NUMNS=30

while getopts :t:p:g:d:o:h opt; do
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
		o)
			OUTFILE=$OPTARG
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

TMPDIR=$(mktemp -d --tmpdir=$TMPOPTION gatkcaller_tmp_XXXXXX)

#Some variables
REFERENCEFILE=$1
FILEIN=$2

###Below uses GATK to do some analysis.
DEDUPLIFIEDBAM=$(mktemp --tmpdir=$TMPDIR --suffix=.bam dedup_XXX)
METRICFILE=$(mktemp --tmpdir=$TMPDIR --suffix=.txt metrics_XXX)
REFERENCEDICT=${REFERENCEFILE%.*}.dict
FULLINTERVALS=$(mktemp --tmpdir=$TMPDIR --suffix=.interval_list fullIntervals_XXX)
SCATTEREDINTERVALDIR=$(mktemp -d --tmpdir=$TMPDIR scatteredIntervals_XXXXXX)
SCATTEREDFIRSTCALLDIR=$(mktemp -d --tmpdir=$TMPDIR scattered_first_calls_XXX)
SUFFIXES=$(seq -f %02.0f 0 $((CORES-1)))
SCATTEREDFIRSTCALLS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i mktemp --tmpdir=$SCATTEREDFIRSTCALLDIR --suffix=.vcf first_calls_{}_XXXXXX)
CMDFIRSTCALLS=$(echo $SCATTEREDFIRSTCALLS | tr ' ' '\n' | xargs -i echo -V {})
JOINEDFIRSTCALLS=$(mktemp --tmpdir=$TMPDIR --suffix=.vcf first_calls_XXX)
RECALIBRATEDBAM=$(mktemp --tmpdir=$TMPDIR --suffix=.bam recal_XXX)
SCATTEREDOUTCALLDIR=$(mktemp -d --tmpdir=$TMPDIR scattered_output_calls_XXX)
SCATTEREDOUTCALLS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i mktemp --tmpdir=$SCATTEREDOUTCALLDIR --suffix=.vcf out_call_{}_XXXXXX)
CMDOUTCALLS=$(echo $SCATTEREDOUTCALLS | tr ' ' '\n' | xargs -i echo -V {})

# REALIGNERINTERVALPREFIX=${TMPDIR}/tmp_intervals_

# INTERVALS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i echo ${REALIGNERINTERVALPREFIX}{}.intervals )
# REALIGNEDBAMS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i echo ${REALIGNERINTERVALPREFIX}{}.bam)
# SORTEDBAMS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i echo ${REALIGNERINTERVALPREFIX}srt_{}.bam)
# REALIGNEDMERGEDBAM=$(mktemp --tmpdir=$TMPDIR --suffix=.bam realigned_XXX)
# RECALDATATABLE=$(mktemp --tmpdir=$TMPDIR --suffix=.table recal_data_XXX)

$PICARD MarkDuplicates INPUT=$FILEIN OUTPUT=$DEDUPLIFIEDBAM METRICS_FILE=$METRICFILE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 || exit 1

if [ ! -e ${REFERENCEFILE}.fai ]; then
	samtools faidx $REFERENCEFILE || exit 1
fi

$PICARD CreateSequenceDictionary REFERENCE=${REFERENCEFILE} OUTPUT=${REFERENCEDICT}
$PICARD BuildBamIndex INPUT=${DEDUPLIFIEDBAM}

#GATK TO RECALIBRATE QUAL SCORES + CALL VARIANTS
$PICARD ScatterIntervalsByNs R=${REFERENCEFILE} OT=ACGT MAX_TO_MERGE=${NUMNS} O=${FULLINTERVALS}
$PICARD IntervalListTools I=${FULLINTERVALS} SCATTER_COUNT=$CORES O=${SCATTEREDINTERVALDIR}
SCATTEREDINTERVALS=$(find ${SCATTEREDINTERVALDIR} -name '*.interval_list')
parallel --halt 2 java -jar ${GATK} -T HaplotypeCaller -R ${REFERENCEFILE} -I $DEDUPLIFIEDBAM -L {1} -stand_call_conf 50 -stand_emit_conf 50 -ploidy 2 -o {2} ::: $SCATTEREDINTERVALS :::+ $SCATTEREDFIRSTCALLS || exit 1
java -cp ${GATK} org.broadinstitute.gatk.tools.CatVariants -R ${REFERENCEFILE} -out ${JOINEDFIRSTCALLS} ${CMDFIRSTCALLS} -assumeSorted
rm $SCATTEREDFIRSTCALLS
java -jar ${GATK} -T BaseRecalibrator -nct $CORES -I $DEDUPLIFIEDBAM -R ${REFERENCEFILE} --knownSites $JOINEDFIRSTCALLS -o $RECALDATATABLE || exit 1
rm $JOINEDFIRSTCALLS
java -jar ${GATK} -T PrintReads -nct $CORES -I $DEDUPLIFIEDBAM -R ${REFERENCEFILE} -BQSR $RECALDATATABLE -EOQ -o $RECALIBRATEDBAM || exit 1
rm $DEDUPLIFIEDBAM $RECALDATATABLE
parallel --halt 2 java -jar ${GATK} -T HaplotypeCaller -R ${REFERENCEFILE} -I $RECALIBRATEDBAM -L {1} -ploidy 2 -o {2} ::: $SCATTEREDINTERVALS :::+ $SCATTEREDOUTCALLS || exit 1
rm $SCATTEREDINTERVALS $RECALIBRATEDBAM
java -cp ${GATK} org.broadinstitute.gatk.tools.CatVariants -R ${REFERENCEFILE} -assumeSorted -out ${OUTFILE} ${CMDOUTCALLS}

rm -rf $TMPDIR

exit 0
