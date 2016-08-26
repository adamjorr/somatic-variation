#!usr/bin/bash
##TODO: REMOVE ALIGNMENT COMPONENTS, UPDATE DOCS + VARIABLES
#aligntoref.sh reference.fasta data.bam
#Takes reference and aligns all .fastq files in any subdirectory and calls SNPs with GATK.
#requires samtools, picard-tools, the GATK, and bowtie2

#http://tldp.org/LDP/abs/html/string-manipulation.html is a great guide for manipulating strings in bash
USAGE="Usage: $0 [-t THREADS] [-p PICARD_CMD] [-g  GATK_PATH] reference.fasta data.bam"

#Here are some things you might want to change:
PICARD="picard" #How do I call picard on this system?
GATK=~/bin/GenomeAnalysisTK.jar #Location of your GATK jar
CORES=48

while getopts :t:p:g:h opt; do
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



#Some variables
REFERENCEFILE=$1
FILEIN=$2

###Below uses GATK to do some analysis.



$PICARD MarkDuplicates INPUT=$FILEIN OUTPUT=dedup.bam METRICS_FILE=metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 || exit 1

if [ ! -e ${REFERENCEFILE}.fai ]; then
	samtools faidx $REFERENCEFILE || exit 1
fi

$PICARD CreateSequenceDictionary REFERENCE=${REFERENCEFILE} OUTPUT=${REFERENCEFILE%.*}.dict
$PICARD BuildBamIndex INPUT=dedup.bam

#Now we use GATK to recalibrate our quality scores and give us a VCF.
java -jar ${GATK} -T RealignerTargetCreator -nt $CORES -R ${REFERENCEFILE} -I dedup.bam -o forIndelAligner.intervals || exit 1
split -d -n l/${CORES} --additional-suffix=.intervals forIndelAligner.intervals tmp_intervals_ || exit 1
parallel java -jar ${GATK} -T IndelRealigner -R ${REFERENCEFILE} -I dedup.bam -targetIntervals tmp_intervals_{}.intervals -o tmp_realigned_{}.bam ::: $(seq -f %02.0f 0 $((CORES-1))) || exit 1
echo tmp_realigned_*.bam  | xargs -d' ' -n 1 -i samtools sort -@ ${CORES} -m 2G -o srt_{} {} || exit 1
samtools merge -@ ${CORES} -c -p realigned.bam srt_tmp_realigned_*.bam || exit 1
java -jar ${GATK} -T UnifiedGenotyper -nt $CORES -I realigned.bam -R ${REFERENCEFILE} -stand_call_conf 50 -stand_emit_conf 50 -ploidy 2 -glm BOTH -o first-calls.vcf || exit 1
java -jar ${GATK} -T BaseRecalibrator -nct $CORES -I realigned.bam -R ${REFERENCEFILE} --knownSites first-calls.vcf -o recal_data.table || exit 1
java -jar ${GATK} -T PrintReads -nct $CORES -I realigned.bam -R ${REFERENCEFILE} -BQSR recal_data.table -EOQ -o recal.bam || exit 1
java -jar ${GATK} -T UnifiedGenotyper -nt $CORES -I recal.bam -R ${REFERENCEFILE} -ploidy 2 -glm BOTH -o var-calls.vcf || exit 1

#Clean up some things
rm aligner* forIndelAligner.intervals realigned.bam first-calls.vcf recal_data.table dedup.bam tmp_intervals_* tmp_realigned_* srt_tmp_realigned_*

exit 0
