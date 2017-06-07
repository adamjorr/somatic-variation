#!/usr/bin/bash
# This is the NCI HPC version which marks duplicates for each BAM file on a single node
##TODO: REMOVE ALIGNMENT COMPONENTS, UPDATE DOCS + VARIABLES
#aligntoref.sh reference.fasta data.bam
#Takes reference and aligns all .fastq files in any subdirectory and calls SNPs with GATK.
#requires samtools, picard-tools, the GATK, and bowtie2

#http://tldp.org/LDP/abs/html/string-manipulation.html is a great guide for manipulating strings in bash
USAGE="Usage: $0 -i bamdir"

#Here are some things you might want to change:
PICARD="java -jar ~/apps/bin/picard.jar" #How do I call picard on this system? on NCI make sure you load the Java module first

while getopts :o:i:h opt; do
        case $opt in
                i)
                        INDIR=$OPTARG
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

if [ $# -eq 0 ]; then                   #if we forget arguments
        echo $USAGE >&2 #remind us
        exit 1                          #and exit with error
fi

if [ "INDIR" == "" ]; then
        echo $USAGE >&2
        echo "Input BAMs directory required." >&2
        exit 1
fi


DEDUPDIR=${INDIR}/dedup
MERGEDDIR=${DEDUPDIR}/merge

mkdir -m777 -p $DEDUPDIR
mkdir -m777 -p $MERGEDDIR

# submit each BAM deduplication job to a cluster node
JOBIDLIST=""
for BAM in $(find $INDIR/*.bam); do
        echo "qsub -v BAM=${BAM},OUTDIR=${DEDUPDIR} picard_dedup.job"
        JOBIDLIST=${JOBIDLIST}:$(qsub -v BAM=${BAM},OUTDIR=$DEDUPDIR HPC_jobs/picard_dedup.job)		# submit a dedup job and record the job id in a string
#        sleep 1
done

# detect that all dedup jobs have completed successfully, then merge all the samples' BAMs into one, then index the big bam
qsub -W depend=afterok${JOBIDLIST} -v DEDUPDIR=${DEDUPDIR},MERGEDDIR=$MERGEDDIR HPC_jobs/bam_merge_index.job

exit 0
