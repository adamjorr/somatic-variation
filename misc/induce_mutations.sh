#!etc/bash

############################################################################################################
# -- written by David Kainer (david.kainer@anu.edu.au)
# This script takes a list of genomic positions (from a reference genome) and fake alternate alleles.
# For a given sample, it then finds all BAM alignments that map across those positions and pulls out the
# relevant reads from that sample's Fastq files (both R1 and R2). It then reqrites the fastq files with
# the new alternate alleles at the correct positions in the relevant reads, taking into account the fact
# that reverse strand reads will be reverse complemented when aligned.
#
# Why do this?
#
# If we place, say, 1000 fake SNPs into the Fastq data and then apply the variant calling pipeline, we can
# evaluate how many of those SNPs are detected and how many are discarded under varying filter settings. This
# gives us a false negative rate for the pipeline.
#############################################################################################################

module load samtools

MUTFILE=/short/xf1/Epauc/test/mutations.txt
REF=/short/xf1/Epauc/pseudoref/RB7_5.fa
BAMDIR=/short/xf1/Epauc/bam/dedup
FQDIR=/short/xf1/Epauc/test
#FQDIR=/short/xf1/Epauc/raw/SN877_merged_paired
/

BAMID=RL41
FQID=RL41_S1
#SNPPOS=1248712
#SNP='='

# handle the R1 reads
sedcommand=""
while read -r CHROM SNPPOS SNP throwaway; 	# for each mutation in the input file, edit the appropriate fastq reads at the right spot!
do
	echo "inducing mutation in R1 Fastq at:" $CHROM $SNPPOS $SNP
	echo "======================================================"

	refbase=$(samtools faidx $REF ${CHROM}:${SNPPOS}-${SNPPOS} | tail -n 1)
	echo "refbase: " $refbase
	# extract the reads which need to be updated for this sample
	# we need to get both forward and reverse reads and treat them differently because reverse have been revcomped in the BAM

	########### reads came from R1 and have been revcomped in the BAM (i.e. from reverse strand)
	READS=$(samtools view -h -f 64 -f 16 $BAMDIR/${BAMID}.dedup.bam ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
 	| awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')

	echo "$READS"

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
 	 	linenum=$(grep -m 1 -nw $readid $FQDIR/${FQID}_R1.fastq | cut -f1 -d':')
 	 	linenum=$((linenum + 1))
 	 	echo "FASTQ line $linenum"

 	 	sedcommand="$sedcommand -e ${linenum}s/${datarc}/${datamod2}/"
	done <<< "$READS"


	######## reads came from R1 and have NOT been revcomped in the BAM (i.e. from forward strand)
	READS=$(samtools view -h -f 64 -F 16 $BAMDIR/${BAMID}.dedup.bam ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
 	| awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')

	echo "$READS"

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

 	 	linenum=$(grep -m 1 -nw $readid $FQDIR/${FQID}_R1.fastq | cut -f1 -d':')
 	 	linenum=$((linenum + 1))
 	 	echo "FASTQ line $linenum"

 	 	sedcommand="$sedcommand -e ${linenum}s/${data}/${datamod}/"
	done <<< "$READS"

done < $MUTFILE

############# execute the actual fastq editing
echo "sed -i ${sedcommand} ${FQDIR}/${FQID}_R1.fastq"
sed -i ${sedcommand} ${FQDIR}/${FQID}_R1.fastq


# handle the R2 reads
sedcommand=""
while read -r CHROM SNPPOS SNP throwaway;
do
        echo "inducing mutation in R2 Fastq at:" $CHROM $SNPPOS $SNP
        echo "======================================================"

	refbase=$(samtools faidx $REF ${CHROM}:${SNPPOS}-${SNPPOS} | tail -n 1)
        echo "refbase: " $refbase

	########### reads came from R2 and have been revcomped in the BAM (i.e. from reverse strand)
	READS=$(samtools view -h -f 128 -f 16 $BAMDIR/${BAMID}.dedup.bam ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
 	| awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')

	echo "$READS"

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

 		linenum=$(grep -m 1 -nw $readid $FQDIR/${FQID}_R2.fastq | cut -f1 -d':')
 		linenum=$((linenum + 1))
 		echo "FASTQ line $linenum"

 		sedcommand="$sedcommand -e ${linenum}s/${datarc}/${datamod2}/"
	done <<< "$READS"


	######## reads came from R2 and have NOT been revcomped in the BAM (i.e. from forward strand)
	READS=$(samtools view -h -f 128 -F 16 $BAMDIR/${BAMID}.dedup.bam ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
 	| awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')

	echo "$READS"
	# NOTE: placing quotes around the variable means that the echo output is AS IS.

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

 		linenum=$(grep -m 1 -nw $readid $FQDIR/${FQID}_R2.fastq | cut -f1 -d':')
 		linenum=$((linenum + 1))
 		echo "FASTQ line $linenum"

 		sedcommand="$sedcommand -e ${linenum}s/${data}/${datamod}/"
	done <<< "$READS"

done < $MUTFILE

############# execute the actual fastq editing
echo "sed -i ${sedcommand} ${FQDIR}/${FQID}_R2.fastq"
sed -i ${sedcommand} ${FQDIR}/${FQID}_R2.fastq

