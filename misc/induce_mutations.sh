\#!/bin/bash

############################################################################################################
# -- written by David Kainer (david.kainer@anu.edu.au)
# This script takes a list of genomic positions (from a reference genome) and fake alternate alleles.
# For a given sample, it then finds all BAM alignments that map across those positions and pulls out the
# relevant reads from that sample's Fastq files (both R1 and R2). It then rewrites the fastq files with
# the new alternate alleles at the correct positions in the relevant reads, taking into account the fact
# that reverse strand reads will be reverse complemented when aligned, and taking into account soft-clipping
#
# Why do this?
#
# If we place, say, 1000 fake SNPs into the Fastq data and then apply the variant calling pipeline, we can
# evaluate how many of those SNPs are detected and how many are discarded under varying filter settings. This
# gives us a false negative rate for the pipeline.
#############################################################################################################

module load samtools

#MUTFILE=/short/xf1/Epauc/test/mutations.txt
#REF=/short/xf1/Epauc/pseudoref/RB7_5.fa
#BAMDIR=/short/xf1/Epauc/bam/dedup
#FQDIR=/short/xf1/Epauc/test
#FQDIR=/short/xf1/Epauc/raw/SN877_merged_paired

#BAMID=RL41
#FQID=RL41_S1
#SNPPOS=1248712
#SNP='='

#########################  get options and their arguments from command line

usage()
{
	echo "-b <bam file>.	The BAM generated for this sample when its original FASTQ was aligned to the reference genome";
	echo "-f <R1 FASTQ>.	The R1 fastq file. i.e. forward reads from illumina PE sequencing";
	echo "-r <R2 FASTQ>.    The R2 fastq file. i.e. reverse reads from illumina PE sequencing";
	echo "-m <mutations>.	A tab delimited file of the mutations to induce in this sample's fastq. Has 3 un-named columns: chromosome, bp position, new allele. e.g. Chr01	228934	C";
	echo "-g <ref genome>.	The reference genome that was used in making the original BAM file"
}

while getopts "b:f:r:m:g:h" opt;
do
	case "${opt}" in
		h ) usage;  exit	;;
		b ) BAM=$OPTARG 	;;
		f ) R1=$OPTARG 		;;
		r ) R2=$OPTARG          ;;
		m ) MUTFILE=$OPTARG     ;;
		g ) REF=$OPTARG		;;
		: ) echo "Option -"$OPTARG" requires an argument" >&2;	 exit 1	 ;;
		* ) echo "invalid option $OPTARG" >&2;   exit    ;;
	esac
done

if [ $OPTIND -lt 10 ]; then echo "all options must be set by the user"; usage; exit; fi
if [ $OPTIND -eq 1 ]; then echo "No options were passed. We really must have options. We simply can't go on without options. Please...options"; exit; fi



# handle the R1 reads
sedcommand=""
while read -r CHROM SNPPOS SNP throwaway; 	# for each mutation in the input file, edit the appropriate fastq reads at the right spot! Getting the right spot is the hard bit :O
do
	echo "inducing mutation in R1 Fastq at:" $CHROM $SNPPOS $SNP
	echo "======================================================"

	refbase=$(samtools faidx $REF ${CHROM}:${SNPPOS}-${SNPPOS} | tail -n 1)
	echo "refbase: " $refbase
	# extract the reads which need to be updated for this sample
	# we need to get both forward and reverse reads and treat them differently because reverse have been revcomped in the BAM

	########### reads came from R1 and have been revcomped in the BAM (i.e. from reverse strand)
#	READS=$(samtools view -h -f 64 -f 16 $BAMDIR/${BAMID}.dedup.bam ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
# 	| awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')

	READS=$(samtools view -h -f 64 -f 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
        | awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')


	echo "$READS"

	if [ -n "$READS" ]
	then
		# for each read to be updated, find it in the fastq and edit the base at the correct offset
		while read -r snppos offset base readid readstart data throwaway; do
 	 		echo "read to edit: " $snppos $offset $base $readid $f5 $data;
 	 		echo "original:" $data

 	 		datarc=$(echo $data | grep '^[ATCG]' - | rev | tr ATCG TAGC) 

 	 		# put the SNP allele into ~50% of the reads and the reference allele into ~50%. This makes it a heterozygous SNP
			if [ $((RANDOM % 2)) -eq 0 ]
                	then
		 	datamod=$(echo $data | sed s/./${SNP}/${offset})                     # replace the specific base at offset with the simulated mutation
 	 		else
		 	datamod=$(echo $data | sed s/./${refbase}/${offset})
			fi

			datamod2=$(echo $datamod | grep '^[ATCG]' - | rev | tr ATCG TAGC)               # reverse complement the alignment
 	 		echo "proposed:" $datamod2

	 		# search fastq file for the read ID and note the line number
		 	# linenum=$(grep -m 1 -nw $readid $FQDIR/${FQID}_R1.fastq | cut -f1 -d':')
			linenum=$(grep -m 1 -nw $readid $R1 | cut -f1 -d':')
 	 		linenum=$((linenum + 1))
 	 		echo "FASTQ line $linenum"

 	 		sedcommand="$sedcommand -e ${linenum}s/${datarc}/${datamod2}/"
		done <<< "$READS"
	fi

	######## reads came from R1 and have NOT been revcomped in the BAM (i.e. from forward strand)
	READS=$(samtools view -h -f 64 -F 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
 	| awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')

	echo "$READS"

	if [ -n "$READS" ]
	then
		while read -r snppos offset base readid readstart data throwaway; do
 	 		echo "read to edit: " $snppos $offset $base $readid $f5 $data;
 	 		# search fastq file of that sample for the read ID

	 		echo "original:" $data

			# put the SNP allele into ~50% of the reads and the reference allele into ~50%. This makes it a heterozygous SNP
                	if [ $((RANDOM % 2)) -eq 0 ]
                	then
                 	datamod=$(echo $data | sed s/./${SNP}/${offset})                     # replace the specific base at offset with the simulated mutation
                	else
                 	datamod=$(echo $data | sed s/./${refbase}/${offset})
                	fi

 	 		echo "proposed:" $datamod

 	 		linenum=$(grep -m 1 -nw $readid $R1 | cut -f1 -d':')
 	 		linenum=$((linenum + 1))
 	 		echo "FASTQ line $linenum"

 	 		sedcommand="$sedcommand -e ${linenum}s/${data}/${datamod}/"
		done <<< "$READS"
	fi
done < $MUTFILE

############# execute the actual fastq editing
filename=$(basename $R1)
fileout="${filename%.*}.mut.fastq"
echo "sed ${sedcommand} ${R1} > $fileout"
sed ${sedcommand} ${R1} > $fileout 
#sed -i ${sedcommand} ${FQDIR}/${FQID}_R1.fastq


# handle the R2 reads
sedcommand=""
while read -r CHROM SNPPOS SNP throwaway;
do
        echo "inducing mutation in R2 Fastq at:" $CHROM $SNPPOS $SNP
        echo "======================================================"

	refbase=$(samtools faidx $REF ${CHROM}:${SNPPOS}-${SNPPOS} | tail -n 1)
        echo "refbase: " $refbase

	########### reads came from R2 and have been revcomped in the BAM (i.e. from reverse strand)
	READS=$(samtools view -h -f 128 -f 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
 	| awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')

	echo "$READS"
	if [ -n "$READS" ]
	then
		while read -r snppos offset base readid readstart data throwaway; do
 			echo "read to edit: " $snppos $offset $base $readid $f5 $data;
 			# search fastq file of that sample for the read ID

 			echo "original:" $data
 			datarc=$(echo $data | grep '^[ATCG]' - | rev | tr ATCG TAGC)

			# put the SNP allele into ~50% of the reads and the reference allele into ~50%. This makes it a heterozygous SNP
                	if [ $((RANDOM % 2)) -eq 0 ]
                	then
                 	datamod=$(echo $data | sed s/./${SNP}/${offset})                     # replace the specific base at offset with the simulated mutation
                	else
                 	datamod=$(echo $data | sed s/./${refbase}/${offset})
                	fi

			datamod2=$(echo $datamod | grep '^[ATCG]' - | rev | tr ATCG TAGC)               # reverse complement the alignment
 			echo "proposed:" $datamod2

 			linenum=$(grep -m 1 -nw $readid $R2 | cut -f1 -d':')
 			linenum=$((linenum + 1))
 			echo "FASTQ line $linenum"

 			sedcommand="$sedcommand -e ${linenum}s/${datarc}/${datamod2}/"
		done <<< "$READS"
	fi

	######## reads came from R2 and have NOT been revcomped in the BAM (i.e. from forward strand)
	READS=$(samtools view -h -f 128 -F 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
 	| awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')

	echo "$READS"
	# NOTE: placing quotes around the variable means that the echo output is AS IS.
	if [ -n "$READS" ]
	then
		while read -r snppos offset base readid readstart data throwaway; do
 			echo "read to edit: " $snppos $offset $base $readid $f5 $data;
 			# search fastq file of that sample for the read ID

 			echo "original:" $data

 			# put the SNP allele into ~50% of the reads and the reference allele into ~50%. This makes it a heterozygous SNP
                	if [ $((RANDOM % 2)) -eq 0 ]
                	then
                 	datamod=$(echo $data | sed s/./${SNP}/${offset})                     # replace the specific base at offset with the simulated mutation
                	else
                 	datamod=$(echo $data | sed s/./${refbase}/${offset})
                	fi

 			echo "proposed:" $datamod

 			linenum=$(grep -m 1 -nw $readid $R2 | cut -f1 -d':')
 			linenum=$((linenum + 1))
 			echo "FASTQ line $linenum"

 			sedcommand="$sedcommand -e ${linenum}s/${data}/${datamod}/"
		done <<< "$READS"
	fi
done < $MUTFILE

############# execute the actual fastq editing
filename=$(basename $R2)
fileout="${filename%.*}.mut.fastq"
echo "sed ${sedcommand} ${R2} > $fileout"
sed ${sedcommand} ${R2} > $fileout

