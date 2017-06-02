#!/usr/bin/bash
# This is the NCI HPC version which marks duplicates for each BAM file on a single node
##TODO: REMOVE ALIGNMENT COMPONENTS, UPDATE DOCS + VARIABLES
#aligntoref.sh reference.fasta data.bam
#Takes reference and aligns all .fastq files in any subdirectory and calls SNPs with GATK.
#requires samtools, picard-tools, the GATK, and bowtie2

#http://tldp.org/LDP/abs/html/string-manipulation.html is a great guide for manipulating strings in bash
USAGE="Usage: $0 [-p PICARD_CMD] [-d TMPDIR] [-o OUTDIR] -i bamdir"

#Here are some things you might want to change:
PICARD="java -jar ~/apps/bin/picard.jar" #How do I call picard on this system? on NCI make sure you load the Java module first

while getopts :p:o:r:i:h opt; do
        case $opt in
                p)
                        PICARD=$OPTARG
                        ;;
                o)
                        OUTDIR=$OPTARG
                        ;;
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


if [ $# -ne 0 ]; then                   #if we forget arguments
        echo $USAGE >&2 #remind us
        exit 1                          #and exit with error
fi

if [ "INDIR" == "" ]; then
        echo $USAGE >&2
        echo "Input BAMs directory required." >&2
        exit 1
fi


mkdir -p $INDIR/dedup
DEDUPDIR=$INDIR/dedup

# submit each BAM deduplication job to a cluster node
for BAM in $(find $INDIR/*.bam); do
        DEDUPLIFIEDBAM=$(mktemp --tmpdir=$TMPDIR --suffix=.bam dedup_XXX)
        echo "qsub -v BAM=${BAM},OUTDIR=${DEDUPDIR} picard_dedup.job"
        qsub -v BAM=${BAM},OUTDIR=${DEDUPDIR} HPC_jobs/picard_dedup.job
        sleep 1
done

exit 0
