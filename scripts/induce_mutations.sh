#!/bin/bash

#########################################################################################################################
# -- written by David Kainer (david.kainer@anu.edu.au)
# This script takes a list of genomic positions (from a reference genome), their original GT and a proposed GT (mutation).
# For a given sample, it then finds in the BAM all alignments that map across those positions and pulls out the
# corresponding reads from that sample's Fastq files (both R1 and R2). It then rewrites the fastq files with
# the new alternate alleles at the correct positions in the relevant reads, taking into account the fact
# that reverse strand reads were reverse complemented when aligned, and taking into account soft-clipping and indels.
#
# Why do this?
#
# If we place, say, 1000 fake SNPs into the Fastq data and then apply a variant calling pipeline, we can
# evaluate how many of those SNPs are detected and how many are discarded under varying filter settings. This
# gives us a false negative rate for the pipeline.
#######################################################################################################################

SCRIPTDIR="$(dirname "$(readlink -f "$0")")"

#########################  get options and their arguments from command line

usage()
{
        echo "-b <bam file>.    The BAM generated for this sample when its original FASTQ was aligned to. Must have .bai index too. FULL path"
        echo "-f <R1 FASTQ>.    The R1 fastq file. i.e. forward reads from illumina PE sequencing. FULL path. NOT .gz"
        echo "-r <R2 FASTQ>.    The R2 fastq file. i.e. reverse reads from illumina PE sequencing. FULL path. NOT .gz"
        echo "-m <mutations>.   A tab delimited file of the mutations to induce in this sample's fastq. FULL path. Expected fields (no headers): CHROM POS REF ALT NEWREF NEWALT"
        echo "-g <ref genome>.  The reference genome that was used in making the original BAM file. Must have .fai index too. FULL path"
        echo "-o <output dir>.  The output directory for the edited fastq files. FULL path"
}


while getopts "b:f:r:m:g:o:h" opt;
do
        case "${opt}" in
                h ) usage;  exit        ;;
                b ) BAM=$OPTARG         ;;
                f ) R1=$OPTARG          ;;
                r ) R2=$OPTARG          ;;
                m ) MUTFILE=$OPTARG     ;;
                g ) GENOME=$OPTARG         ;;
                o ) OUTDIR=$OPTARG      ;;
                : ) echo "Option -"$OPTARG" requires an argument" >&2;   exit 1  ;;
                * ) echo "invalid option $OPTARG" >&2;   exit    ;;
        esac
done

if [ $OPTIND -lt 12 ]; then echo "all options MUST be set by the user"; usage; exit; fi
if [ $OPTIND -eq 1 ]; then echo "No options were passed. We really must have options. We simply can't go on"; usage; exit; fi

# handle the R1 reads
count=0
sedcommand=""

while read -r CHROM SNPPOS REF ALT NEWREF NEWALT throwaway;     # for each mutation in the input file
do
        count=$((count+1))
        echo -e "\n $count - inducing mutation $NEWREF $NEWALT in R1 Fastq at:" $CHROM $SNPPOS >> $OUTDIR/$(basename $MUTFILE).log
        echo "===========================================================================" >> $OUTDIR/$(basename $MUTFILE).log

        refbase=$(samtools faidx $GENOME ${CHROM}:${SNPPOS}-${SNPPOS} | tail -n 1)
        echo "ref genome base: " $refbase "| ROOT GT: " $REF $ALT  >> $OUTDIR/$(basename $MUTFILE).log

	# extract the reads which need to be updated for this sample
        # we need to get both forward and reverse reads and treat them differently because reverse have been revcomped in the BAM

        ########### reads came from R1 and have been revcomped in the BAM (i.e. from reverse strand)
	READS=$(samtools view -h -f 64 -f 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
	| awk -v sdir=$SCRIPTDIR -v pos=$SNPPOS 'BEGIN {OFS=FS="\t"} ; {n=split($10,a,"")}; {cigar=$6; offset=(pos-$4)+1; cmd="perl -s "sdir"/cigarpos.pl -p="offset" -c="cigar; cmd | getline adjOffset; close(cmd);  print pos, offset, adjOffset, a[adjOffset], $1, $4, cigar, $10}')

	#echo "# of R1 reverse strand reads to edit: " $(echo -n "$READS" | grep -c '^')  >> $OUTDIR/logfile.txt

        if [ -n "$READS" ]
        then
                # for each read to be updated, find it in the fastq and edit the base at the correct offset
		logstring=""
                while read -r snppos offset adjOffset base readid readstart cigar data throwaway;
                do
			# Check case where perl script fails to provide an adjusted offset because the  unadjusted offset lands in a deletion. So skip the read in that case!
                        [ "$adjOffset" -eq -1 ] && continue

                        logstring="$logstring \n read to edit:  $snppos $offset $adjOffset $base $readid $readstart $cigar"
                        logstring="$logstring \n alignment:\t $data"

			datamod=$data
			# if site was homozygous and the read's base is NOT an error, then mutate it with a 50% chance to make the site heterozygous
			if [ "$REF" = "$ALT" ] && [ "$base" = "$ALT" ]; then
				if [ $((RANDOM % 2)) -eq 0 ]; then
				datamod=$(echo $data | sed s/./${NEWALT}/${adjOffset})
				fi
			# but if the site was already heterozygous...
			elif [ "$REF" != "$ALT" ] && [ "$ALT" != "$NEWALT" ] && [ "$base" = "$ALT" ]; then	#...and if the ALT allele is to be mutated, and if the read's base is an ALT base, then mutate it!
				datamod=$(echo $data | sed s/./${NEWALT}/${adjOffset})
			elif [ "$REF" != "$ALT" ] && [ "$REF" != "$NEWREF" ] && [ "$base" = "$REF" ]; then	#...or if the REF allele is to be mutated,  and if the read's base is a  REF base, then mutate it!
				datamod=$(echo $data | sed s/./${NEWREF}/${adjOffset})
			fi

                        logstring="$logstring \n proposed:\t  $datamod"

			## ...only create a sed command if an actual change to this read is occuring. i.e. if $data != $datamod
			if [ "$datamod" != "$data" ]
			then
                        	datamod2=$(echo $datamod | grep '^[ATCGN]' - | rev | tr ATCG TAGC)                # reverse complement the modified alignment
				datarc=$(echo $data | grep '^[ATCGN]' - | rev | tr ATCG TAGC)
                        	# search fastq file for the read ID and note the line number
                        	linenum=$(grep -m 1 -nw $readid $R1 | cut -f1 -d':')
                        	linenum=$((linenum + 1))
                        	#echo "FASTQ line $linenum"  >> $OUTDIR/logfile.txt
                        	logstring="$logstring \n sed command:\t ${linenum}s/${datarc}/${datamod2}/"
                        	sedcommand="$sedcommand ${linenum}s/${datarc}/${datamod2}/\n"
			fi
                done <<< "$READS"
		echo -e $logstring >> $OUTDIR/$(basename $MUTFILE).log
        fi

	######## reads came from R1 and have NOT been revcomped in the BAM (i.e. from forward strand)
        READS=$(samtools view -h -f 64 -F 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
        | awk -v sdir=$SCRIPTDIR -v pos=$SNPPOS 'BEGIN {OFS=FS="\t"} ; {n=split($10,a,"")}; {cigar=$6; offset=(pos-$4)+1; cmd="perl -s "sdir"/cigarpos.pl -p="offset" -c="cigar; cmd | getline adjOffset; close(cmd);  print pos, offset, adjOffset, a[adjOffset], $1, $4, cigar, $10}')

        #output: genomepos  unadjOffset  adjOffset base readid  readstart CIGAR read

        #echo "# of R1 forward strand reads to edit: " $(echo -n "$READS" | grep -c '^')  >> $OUTDIR/logfile.txt

        if [ -n "$READS" ]
        then
		logstring=""
                while read -r snppos offset adjOffset base readid readstart cigar data throwaway;
		do
			[ "$adjOffset" -eq -1 ] && continue        # Check case where perl script fails to provide an adjusted offset when the unadjusted offset it lands in a deletion. So skip the read in that case!

			logstring="$logstring \n read to edit:  $snppos $offset $adjOffset $base $readid $readstart $cigar"
                        logstring="$logstring \n alignment:\t $data"

			datamod=$data
                        # if site was homozygous and the read's base is NOT an error, then mutate it with a 50% chance to make the site heterozygous
                        if [ "$REF" = "$ALT" ] && [ "$base" = "$ALT" ]; then
                                if [ $((RANDOM % 2)) -eq 0 ]; then
                                datamod=$(echo $data | sed s/./${NEWALT}/${adjOffset})
                                fi
                        # but if the site was already heterozygous...
                        elif [ "$REF" != "$ALT" ] && [ "$ALT" != "$NEWALT" ] && [ "$base" = "$ALT" ]; then      #...and if the ALT allele is to be mutated, and if the read's base is an ALT base, then mutate it!
                                datamod=$(echo $data | sed s/./${NEWALT}/${adjOffset})
                        elif [ "$REF" != "$ALT" ] && [ "$REF" != "$NEWREF" ] && [ "$base" = "$REF" ]; then      #...or if the REF allele is to be mutated,  and if the read's base is a  REF base, then mutate it!
                                datamod=$(echo $data | sed s/./${NEWREF}/${adjOffset})
                        fi

                        logstring="$logstring \n proposed:\t $datamod"

			## ...only create a sed command if an actual change to this read is occuring. i.e. if $data != $datamod
                        if [ "$datamod" != "$data" ]
                        then
	                        linenum=$(grep -m 1 -nw $readid $R1 | cut -f1 -d':')
        	                linenum=$((linenum + 1))
                	        #echo "FASTQ line $linenum" >> $OUTDIR/logfile.txt
				logstring="$logstring \n sed command:\t ${linenum}s/${data}/${datamod}/"
                        	sedcommand="$sedcommand ${linenum}s/${data}/${datamod}/\n"
			fi
                done <<< "$READS"
		echo -e $logstring >> $OUTDIR/$(basename $MUTFILE).log
        fi

done < $MUTFILE


############# output the sed command for this chunk of mutations
filename=$(basename $MUTFILE)
fileout="${filename}.R1.sed"
echo -e "${sedcommand}" > ${OUTDIR}/$fileout



# handle the R2 reads
count=0
sedcommand=""
while read -r  CHROM SNPPOS REF ALT NEWREF NEWALT throwaway;
do
        count=$((count+1))
        echo -e "\n $MUTFILE $count - inducing mutation $SNP in R2 Fastq at:" $CHROM $SNPPOS >> $OUTDIR/$(basename $MUTFILE).log
        echo "================================================================" >> $OUTDIR/$(basename $MUTFILE).log

	refbase=$(samtools faidx $GENOME ${CHROM}:${SNPPOS}-${SNPPOS} | tail -n 1)
        echo "$MUTFILE ref genome base: " $refbase "| ROOT GT: " $REF $ALT >> $OUTDIR/$(basename $MUTFILE).log

        ########### reads came from R2 and have been revcomped in the BAM (i.e. from reverse strand)
        READS=$(samtools view -h -f 128 -f 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
        | awk -v sdir=$SCRIPTDIR -v pos=$SNPPOS 'BEGIN {OFS=FS="\t"} ; {n=split($10,a,"")}; {cigar=$6; offset=(pos-$4)+1; cmd="perl -s "sdir"/cigarpos.pl -p="offset" -c="cigar; cmd | getline adjOffset; close(cmd);  print pos, offset,adjOffset, a[adjOffset], $1, $4, cigar, $10}')
	#output: genomepos  unadjOffset  adjOffset base readid  readstart CIGAR read

        #echo "# of R2 reverse strand reads to edit: " $(echo -n "$READS" | grep -c '^') >> $OUTDIR/logfile.txt

        if [ -n "$READS" ]
        then
		logstring=""
                 while read -r snppos offset adjOffset base readid readstart cigar data throwaway;
                        do
                        [ "$adjOffset" -eq -1 ] && continue     # Check case where perl script fails to provide an adjusted offset when the unadjusted offset it lands in a deletion. So skip the read in that case!

			logstring="$logstring \n read to edit:  $snppos $offset $adjOffset $base $readid $readstart $cigar"
                        logstring="$logstring \n alignment:\t $data"

			datamod=$data
                        # if site was homozygous and the read's base is NOT an error, then mutate it with a 50% chance to make the site heterozygous
                        if [ "$REF" = "$ALT" ] && [ "$base" = "$ALT" ]; then
                                if [ $((RANDOM % 2)) -eq 0 ]; then
                                datamod=$(echo $data | sed s/./${NEWALT}/${adjOffset})
                                fi
                        # but if the site was already heterozygous...
                        elif [ "$REF" != "$ALT" ] && [ "$ALT" != "$NEWALT" ] && [ "$base" = "$ALT" ]; then      #...and if the ALT allele is to be mutated, and if the read's base is an ALT base, then mutate it!
                                datamod=$(echo $data | sed s/./${NEWALT}/${adjOffset})
                        elif [ "$REF" != "$ALT" ] && [ "$REF" != "$NEWREF" ] && [ "$base" = "$REF" ]; then      #...or if the REF allele is to be mutated,  and if the read's base is a  REF base, then mutate it!
                                datamod=$(echo $data | sed s/./${NEWREF}/${adjOffset})
                        fi

                        logstring="$logstring \n proposed:\t $datamod"

 			## ...only create a sed command if an actual change to this read is occuring. i.e. if $data != $datamod
                        if [ "$datamod" != "$data" ]
			then
	                        datamod2=$(echo $datamod | grep '^[ATCGN]' - | rev | tr ATCG TAGC)               # reverse complement the modified alignment
				datarc=$(echo $data | grep '^[ATCGN]' - | rev | tr ATCG TAGC)
	                        linenum=$(grep -m 1 -nw $readid $R2 | cut -f1 -d':')
        	                linenum=$((linenum + 1))
                	        #echo "FASTQ line $linenum" >> $OUTDIR/logfile.txt
				logstring="$logstring \n sed command:\t ${linenum}s/${datarc}/${datamod2}/"
                        	sedcommand="$sedcommand ${linenum}s/${datarc}/${datamod2}/\n"
			fi
                done <<< "$READS"
		echo -e $logstring >> $OUTDIR/$(basename $MUTFILE).log
        fi

	######## reads came from R2 and have NOT been revcomped in the BAM (i.e. from forward strand)
        READS=$(samtools view -h -f 128 -F 16 ${BAM} ${CHROM}:${SNPPOS}-${SNPPOS} | grep -v "^@" \
        | awk -v sdir=$SCRIPTDIR -v pos=$SNPPOS 'BEGIN {OFS=FS="\t"} ; {n=split($10,a,"")}; {cigar=$6; offset=(pos-$4)+1; cmd="perl -s "sdir"/cigarpos.pl -p="offset" -c="cigar; cmd | getline adjOffset; close(cmd);  print pos, offset, adjOffset, a[adjOffset], $1, $4, cigar, $10}')
        #output: genomepos  unadjOffset  adjOffset base readid  readstart CIGAR read

        #echo "# of R2 forward strand reads to edit: " $(echo -n "$READS" | grep -c '^') >> $OUTDIR/logfile.txt

        # NOTE: placing quotes around the variable means that the echo output is AS IS.
        if [ -n "$READS" ]
        then
                while read -r snppos offset adjOffset base readid readstart cigar data throwaway;
                do
                        [ "$adjOffset" -eq -1 ] && continue       # Check case where perl script fails to provide an adjusted offset when the unadjusted offset it lands in a deletion. So skip the read in that case!

			logstring="$logstring \n read to edit:  $snppos $offset $adjOffset $base $readid $readstart $cigar"
                        logstring="$logstring \n alignment:\t $data"

			datamod=$data
                        # if site was homozygous and the read's base is NOT an error, then mutate it with a 50% chance to make the site heterozygous
                        if [ "$REF" = "$ALT" ] && [ "$base" = "$ALT" ]; then
                                if [ $((RANDOM % 2)) -eq 0 ]; then
                                datamod=$(echo $data | sed s/./${NEWALT}/${adjOffset})
                                fi
                        # but if the site was already heterozygous...
                        elif [ "$REF" != "$ALT" ] && [ "$ALT" != "$NEWALT" ] && [ "$base" = "$ALT" ]; then      #...and if the ALT allele is to be mutated, and if the read's base is an ALT base, then mutate it!
                                datamod=$(echo $data | sed s/./${NEWALT}/${adjOffset})
                        elif [ "$REF" != "$ALT" ] && [ "$REF" != "$NEWREF" ] && [ "$base" = "$REF" ]; then      #...or if the REF allele is to be mutated,  and if the read's base is a  REF base, then mutate it!
                                datamod=$(echo $data | sed s/./${NEWREF}/${adjOffset})
                        fi

			logstring="$logstring \n proposed:\t $datamod"

 			## ...only create a sed command if an actual change to this read is occuring. i.e. if $data != $datamod
                        if [ "$datamod" != "$data" ]
                        then
	                        linenum=$(grep -m 1 -nw $readid $R2 | cut -f1 -d':')
        	                linenum=$((linenum + 1))
                	        #echo "FASTQ line $linenum" >> $OUTDIR/logfile.txt
				logstring="$logstring \n sed command:\t ${linenum}s/${data}/${datamod}/"
                        	sedcommand="$sedcommand ${linenum}s/${data}/${datamod}/\n"
			fi
                done <<< "$READS"
		echo -e $logstring >> $OUTDIR/$(basename $MUTFILE).log
        fi
done < $MUTFILE

############# output the sed command for this chunk of mutations
filename=$(basename $MUTFILE)
fileout="${filename}.R2.sed"
echo -e "${sedcommand}" > ${OUTDIR}/$fileout
