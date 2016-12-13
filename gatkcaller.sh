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
REALIGNERINTERVALS=$(mktemp --tmpdir=$TMPDIR --suffix=.intervals forIndelAligner_XXX)
REALIGNERINTERVALPREFIX=${TMPDIR}/tmp_intervals_
SUFFIXES=$(seq -f %02.0f 0 $((CORES-1)))
INTERVALS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i echo ${REALIGNERINTERVALPREFIX}{}.intervals )
REALIGNEDBAMS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i echo ${REALIGNERINTERVALPREFIX}{}.bam)
SORTEDBAMS=$(echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i echo ${REALIGNERINTERVALPREFIX}srt_{}.bam)
REALIGNEDMERGEDBAM=$(mktemp --tmpdir=$TMPDIR --suffix=.bam realigned_XXX)
RECALIBRATEDBAM=$(mktemp --tmpdir=$TMPDIR --suffix=.bam recal_XXX)
FIRSTCALLS=$(mktemp --tmpdir=$TMPDIR --suffix=.vcf first_calls_XXX)
RECALDATATABLE=$(mktemp --tmpdir=$TMPDIR --suffix=.table recal_data_XXX)

$PICARD MarkDuplicates INPUT=$FILEIN OUTPUT=$DEDUPLIFIEDBAM METRICS_FILE=$METRICFILE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 || exit 1

if [ ! -e ${REFERENCEFILE}.fai ]; then
	samtools faidx $REFERENCEFILE || exit 1
fi

$PICARD CreateSequenceDictionary REFERENCE=${REFERENCEFILE} OUTPUT=${REFERENCEDICT}
$PICARD BuildBamIndex INPUT=${DEDUPLIFIEDBAM}

#Now we use GATK to recalibrate our quality scores and give us a VCF.
java -jar ${GATK} -T RealignerTargetCreator -nt $CORES -R ${REFERENCEFILE} -I $DEDUPLIFIEDBAM -o $REALIGNERINTERVALS || exit 1
split -d -n l/${CORES} --additional-suffix=.intervals $REALIGNERINTERVALS ${REALIGNERINTERVALPREFIX} || exit 1
parallel --halt 2 java -jar ${GATK} -T IndelRealigner -R ${REFERENCEFILE} -I $DEDUPLIFIEDBAM -targetIntervals ${REALIGNERINTERVALPREFIX}{}.intervals -o ${REALIGNERINTERVALPREFIX}{}.bam ::: $SUFFIXES || exit 1
rm $DEDUPLIFIEDBAM $METRICFILE $REALIGNERINTERVALS
echo $SUFFIXES | tr ' ' '\n' | xargs -n 1 -i samtools sort -@ ${CORES} -T ${TMPDIR}/samtmp -m 2G -o ${REALIGNERINTERVALPREFIX}srt_{}.bam ${REALIGNERINTERVALPREFIX}{}.bam || exit 1
rm $INTERVALS $REALIGNEDBAMS
samtools merge -@ ${CORES} -f -c -p $REALIGNEDMERGEDBAM $SORTEDBAMS || exit 1
samtools index $REALIGNEDMERGEDBAM || exit 1
rm $SORTEDBAMS
java -jar ${GATK} -T UnifiedGenotyper -nt $CORES -I $REALIGNEDMERGEDBAM -R ${REFERENCEFILE} -stand_call_conf 50 -stand_emit_conf 50 -ploidy 2 -glm BOTH -o $FIRSTCALLS || exit 1
java -jar ${GATK} -T BaseRecalibrator -nct $CORES -I $REALIGNEDMERGEDBAM -R ${REFERENCEFILE} --knownSites $FIRSTCALLS -o $RECALDATATABLE || exit 1
rm $FIRSTCALLS
java -jar ${GATK} -T PrintReads -nct $CORES -I $REALIGNEDMERGEDBAM -R ${REFERENCEFILE} -BQSR $RECALDATATABLE -EOQ -o $RECALIBRATEDBAM || exit 1
rm $REALIGNEDMERGEDBAM $RECALDATATABLE
java -jar ${GATK} -T UnifiedGenotyper -nt $CORES -I $RECALIBRATEDBAM -R ${REFERENCEFILE} -ploidy 2 -glm BOTH -o $OUTFILE || exit 1

rm $RECALIBRATEDBAM
rm -rf $TMPDIR

exit 0
