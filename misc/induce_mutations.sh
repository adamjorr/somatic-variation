#!/bin/bash

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

SCRIPTDIR=/short/xf1/src_big/somatic-variation/misc
#MUTFILE=/short/xf1/Epauc/test/mutations.txt
#REF=/short/xf1/Epauc/pseudoref/RB7_5.fa
BAMDIR=/short/xf1/Epauc/bam/dedup
#FQDIR=/short/xf1/Epauc/test
#FQDIR=/short/xf1/Epauc/raw/SN877_merged_paired

#BAMID=RL41
#FQID=RL41_S1
#SNPPOS=1248712
#SNP='='

#########################  get options and their arguments from command line

usage()
{
	echo "-b <bam file>.	The BAM generated for this sample when its original FASTQ was aligned to the reference genome. FULL path";
	echo "-f <R1 FASTQ>.	The R1 fastq file. i.e. forward reads from illumina PE sequencing. FULL path";
	echo "-r <R2 FASTQ>.    The R2 fastq file. i.e. reverse reads from illumina PE sequencing. FULL path";
	echo "-m <mutations>.	A tab delimited file of the mutations to induce in this sample's fastq. Has 3 un-named columns: chromosome, bp position, new allele. e.g. Chr01	228934	C";
	echo "-g <ref genome>.	The reference genome that was used in making the original BAM file"
	echo "-o <output dir>.	The output directory for the edited fastq files"
}

while getopts "b:f:r:m:g:o:h" opt;
do
	case "${opt}" in
		h ) usage;  exit	;;
		b ) BAM=$OPTARG 	;;
		f ) R1=$OPTARG 		;;
		r ) R2=$OPTARG          ;;
		m ) MUTFILE=$OPTARG     ;;
		g ) REF=$OPTARG		;;
		o ) OUTDIR=$OPTARG	;;
		: ) echo "Option -"$OPTARG" requires an argument" >&2;	 exit 1	 ;;
		* ) echo "invalid option $OPTARG" >&2;   exit    ;;
	esac
done

if [ $OPTIND -lt 12 ]; then echo "all options MUST be set by the user"; usage; exit; fi
if [ $OPTIND -eq 1 ]; then echo "No options were passed. We really must have options. We simply can't go on without options. Please...options"; exit; fi



# handle the R1 reads
count=0
sedcommand=""
while read -r CHROM SNPPOS SNP throwaway; 	# for each mutation in the input file, edit the appropriate fastq reads at the right spot! Getting the right spot is the hard bit :O
do
	count=$((count+1))
	echo -e "\n $count - inducing mutation $SNP in R1 Fastq at:" $CHROM $SNPPOS
	echo "========================================================"

	refbase=$(samtools faidx $REF ${CHROM}:${SNPPOS}-${SNPPOS} | tail -n 1)
	echo "refbase: " $refbase
	# extract the reads which need to be updated for this sample
	# we need to get both forward and reverse reads and treat them differently because reverse have been revcomped in the BAM

	########### reads came from R1 and have been revcomped in the BAM (i.e. from reverse strand)
	READS=$(samtools view -h -f 64 -f 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
        | awk -v sdir=$SCRIPTDIR -v pos=$SNPPOS 'BEGIN {OFS=FS="\t"} ; {n=split($10,a,"")}; {cigar=$6; offset=(pos-$4)+1; cmd="perl -s "sdir"/cigarpos.pl -p="offset" -c="cigar; cmd | getline adjOffset; close(cmd);  print pos, offset, adjOffset, a[adjOffset], $1, $4, cigar, $10}')
	#output: genomepos  unadjOffset  adjOffset base readid  readstart CIGAR read


	echo "# of R1 reverse strand reads to edit: " $(echo -n "$READS" | grep -c '^')

	if [ -n "$READS" ]
	then
		# for each read to be updated, find it in the fastq and edit the base at the correct offset
		while read -r snppos offset adjOffset base readid readstart cigar data throwaway;
		do
 	 		# Check for case where perl script gives adjusted offset = -1 (i.e. when the unadjusted offset lands in a deletion). Skip the read in that case!
			[ "$adjOffset" -eq -1 ] && continue

			echo "read to edit: " $snppos $offset $adjOffset $base $readid $readstart $cigar;
 	 		echo "alignment: " $data

			datarc=$(echo $data | grep '^[ATCGN]' - | rev | tr ATCG TAGC)

 	 		# put the SNP allele into ~50% of the reads and the reference allele into ~50%. This makes it a heterozygous SNP
			if [ $((RANDOM % 2)) -eq 0 ]
                	then
		 	datamod=$(echo $data | sed s/./${SNP}/${adjOffset})                     # replace the specific base at offset with the simulated mutation
 	 		else
		 	datamod=$(echo $data | sed s/./${refbase}/${adjOffset})
			fi

			echo "proposed:  " $datamod
			datamod2=$(echo $datamod | grep '^[ATCGN]' - | rev | tr ATCG TAGC)               # reverse complement the alignment

	 		# search fastq file for the read ID and note the line number
		 	# linenum=$(grep -m 1 -nw $readid $FQDIR/${FQID}_R1.fastq | cut -f1 -d':')
			linenum=$(grep -m 1 -nw $readid $R1 | cut -f1 -d':')
 	 		linenum=$((linenum + 1))
 	 		echo "FASTQ line $linenum"
			echo "sed command: -e ${linenum}s/${datarc}/${datamod2}/"
 	 		sedcommand="$sedcommand -e ${linenum}s/${datarc}/${datamod2}/"
		done <<< "$READS"
	fi

	######## reads came from R1 and have NOT been revcomped in the BAM (i.e. from forward strand)
	READS=$(samtools view -h -f 64 -F 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
        | awk -v sdir=$SCRIPTDIR -v pos=$SNPPOS 'BEGIN {OFS=FS="\t"} ; {n=split($10,a,"")}; {cigar=$6; offset=(pos-$4)+1; cmd="perl -s "sdir"/cigarpos.pl -p="offset" -c="cigar; cmd | getline adjOffset; close(cmd);  print pos, offset, adjOffset, a[adjOffset], $1, $4, cigar, $10}')

        #output: genomepos  unadjOffset  adjOffset base readid  readstart CIGAR read

	echo "# of R1 forward strand reads to edit: " $(echo -n "$READS" | grep -c '^')

	if [ -n "$READS" ]
	then
		while read -r snppos offset adjOffset base readid readstart cigar data throwaway;
			do
			[ "$adjOffset" -eq -1 ] && continue        # Check case where perl script fails to provide an adjusted offset when the unadjusted offset it lands in a deletion. So skip the read in that case!

			echo "read to edit: " $snppos $offset $adjOffset $base $readid $readstart $cigar;
	 		echo "alignment: " $data

			# put the SNP allele into ~50% of the reads and the reference allele into ~50%. This makes it a heterozygous SNP
                	if [ $((RANDOM % 2)) -eq 0 ]
                	then
                 	datamod=$(echo $data | sed s/./${SNP}/${adjOffset})                     # replace the specific base at offset with the simulated mutation
                	else
                 	datamod=$(echo $data | sed s/./${refbase}/${adjOffset})
                	fi

 	 		echo "proposed:  " $datamod

 	 		linenum=$(grep -m 1 -nw $readid $R1 | cut -f1 -d':')
 	 		linenum=$((linenum + 1))
 	 		echo "FASTQ line $linenum"
			echo "sed command: -e ${linenum}s/${data}/${datamod}/"
 	 		sedcommand="$sedcommand -e ${linenum}s/${data}/${datamod}/"
		done <<< "$READS"
	fi
done < $MUTFILE

############# execute the actual fastq editing
filename=$(basename $R1)
fileout="${filename%.*}.mut.fastq"
echo "sed ${sedcommand} ${R1} > ${OUTDIR}/$fileout"
sed ${sedcommand} ${R1} > ${OUTDIR}/$fileout 


# handle the R2 reads
count=0
sedcommand=""
while read -r CHROM SNPPOS SNP throwaway;
do
	count=$((count+1))
        echo -e "\n $count - inducing mutation $SNP in R2 Fastq at:" $CHROM $SNPPOS
        echo "========================================================"

	refbase=$(samtools faidx $REF ${CHROM}:${SNPPOS}-${SNPPOS} | tail -n 1)
        echo "refbase: " $refbase

	########### reads came from R2 and have been revcomped in the BAM (i.e. from reverse strand)
	READS=$(samtools view -h -f 128 -f 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
	| awk -v sdir=$SCRIPTDIR -v pos=$SNPPOS 'BEGIN {OFS=FS="\t"} ; {n=split($10,a,"")}; {cigar=$6; offset=(pos-$4)+1; cmd="perl -s "sdir"/cigarpos.pl -p="offset" -c="cigar; cmd | getline adjOffset; close(cmd);  print pos, offset, adjOffset, a[adjOffset], $1, $4, cigar, $10}')

        #output: genomepos  unadjOffset  adjOffset base readid  readstart CIGAR read

	echo "# of R2 reverse strand reads to edit: " $(echo -n "$READS" | grep -c '^')

	if [ -n "$READS" ]
	then
		 while read -r snppos offset adjOffset base readid readstart cigar data throwaway;
			do
			[ "$adjOffset" -eq -1 ] && continue     # Check case where perl script fails to provide an adjusted offset when the unadjusted offset it lands in a deletion. So skip the read in that case!
                        echo "read to edit: " $snppos $offset $adjOffset $base $readid $readstart $cigar;
 			echo "alignment: " $data
 			datarc=$(echo $data | grep '^[ATCGN]' - | rev | tr ATCG TAGC)

			# put the SNP allele into ~50% of the reads and the reference allele into ~50%. This makes it a heterozygous SNP
                	if [ $((RANDOM % 2)) -eq 0 ]
                	then
                 	datamod=$(echo $data | sed s/./${SNP}/${adjOffset})                     # replace the specific base at offset with the simulated mutation
                	else
                 	datamod=$(echo $data | sed s/./${refbase}/${adjOffset})
                	fi

			echo "proposed:  " $datamod
			datamod2=$(echo $datamod | grep '^[ATCGN]' - | rev | tr ATCG TAGC)               # reverse complement the modified alignment

 			linenum=$(grep -m 1 -nw $readid $R2 | cut -f1 -d':')
 			linenum=$((linenum + 1))
 			echo "FASTQ line $linenum"
			echo "sed command: -e ${linenum}s/${datarc}/${datamod2}/"
 			sedcommand="$sedcommand -e ${linenum}s/${datarc}/${datamod2}/"
		done <<< "$READS"
	fi

	######## reads came from R2 and have NOT been revcomped in the BAM (i.e. from forward strand)
	READS=$(samtools view -h -f 128 -F 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
	| awk -v sdir=$SCRIPTDIR -v pos=$SNPPOS 'BEGIN {OFS=FS="\t"} ; {n=split($10,a,"")}; {cigar=$6; offset=(pos-$4)+1; cmd="perl -s "sdir"/cigarpos.pl -p="offset" -c="cigar; cmd | getline adjOffset; close(cmd);  print pos, offset, adjOffset, a[adjOffset], $1, $4, cigar, $10}')
        #output: genomepos  unadjOffset  adjOffset base readid  readstart CIGAR read

	echo "# of R2 forward strand reads to edit: " $(echo -n "$READS" | grep -c '^')

	# NOTE: placing quotes around the variable means that the echo output is AS IS.
	if [ -n "$READS" ]
	then
 		while read -r snppos offset adjOffset base readid readstart cigar data throwaway;
		do
			[ "$adjOffset" -eq -1 ] && continue       # Check case where perl script fails to provide an adjusted offset when the unadjusted offset it lands in a deletion. So skip the read in that case!
                        echo "read to edit: " $snppos $offset $adjOffset $base $readid $readstart $cigar;
 			echo "alignment: " $data

 			# put the SNP allele into ~50% of the reads and the reference allele into ~50%. This makes it a heterozygous SNP
                	if [ $((RANDOM % 2)) -eq 0 ]
                	then
                 	datamod=$(echo $data | sed s/./${SNP}/${adjOffset})                     # replace the specific base at offset with the simulated mutation
                	else
                 	datamod=$(echo $data | sed s/./${refbase}/${adjOffset})
                	fi

 			echo "proposed:  " $datamod

 			linenum=$(grep -m 1 -nw $readid $R2 | cut -f1 -d':')
 			linenum=$((linenum + 1))
 			echo "FASTQ line $linenum"
			echo "sed command: -e ${linenum}s/${data}/${datamod}/"
 			sedcommand="$sedcommand -e ${linenum}s/${data}/${datamod}/"
		done <<< "$READS"
	fi
done < $MUTFILE

############# execute the actual fastq editing
filename=$(basename $R2)
fileout="${filename%.*}.mut.fastq"
echo "sed ${sedcommand} ${R2} > ${OUTDIR}/$fileout"
sed ${sedcommand} ${R2} > ${OUTDIR}/$fileout

