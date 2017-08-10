#!etc/bash

module load samtools

BAMDIR=/short/xf1/Epauc/bam/dedup
FQDIR=/short/xf1/Epauc/test
#FQDIR=/short/xf1/Epauc/raw/SN877_merged_paired


BAMID=RL41
FQID=RL41_S1
SNPPOS=1248712
SNP='='


# first we extract the reads which need to be updated for this sample
# we need to get both forward and reverse reads and treat them differently because reverse have been revcomped in the BAM

sedcommand=""
########### reads came from R1 and have been revcomped in the BAM (i.e. from reverse strand)
READS=$(samtools view -h -f 64 -f 16 $BAMDIR/${BAMID}.dedup.bam Chr01:${SNPPOS}-${SNPPOS} | grep -v "^@" \
 | awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')

echo "$READS"

# for each read to be updated, find it in the fastq and edit the base at the correct offset
while read -r snppos offset base readid readstart data throwaway; do
 echo "read to edit: " $snppos $offset $base $readid $f5 $data;
 # search fastq file of that sample for the read ID

 echo "original:" $data
 datarc=$(echo $data | grep '^[ATCG]' - | rev | tr ATCG TAGC) 
 datamod=$(echo $data | sed s/./${SNP}/${offset})                               # replace the specific base at offset with the simulated mutation
 datamod2=$(echo $datamod | grep '^[ATCG]' - | rev | tr ATCG TAGC)               # reverse complement the alignment
 echo "proposed:" $datamod2

 linenum=$(grep -m 1 -nw $readid $FQDIR/${FQID}_R1.fastq | cut -f1 -d':')
 linenum=$((linenum + 1))
 echo "FASTQ line $linenum"

 sedcommand="$sedcommand -e ${linenum}s/${datarc}/${datamod2}/"
done <<< "$READS"


######## reads came from R1 and have NOT been revcomped in the BAM (i.e. from forward strand)
READS=$(samtools view -h -f 64 -F 16 $BAMDIR/${BAMID}.dedup.bam Chr01:${SNPPOS}-${SNPPOS} | grep -v "^@" \
 | awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')


echo "$READS"

while read -r snppos offset base readid readstart data throwaway; do
 echo "read to edit: " $snppos $offset $base $readid $f5 $data;
 # search fastq file of that sample for the read ID

 echo "original:" $data
 datamod=$(echo $data | sed s/./${SNP}/${offset}) # replace the specific base at offset with the simulated mutation
 echo "proposed:" $datamod

 linenum=$(grep -m 1 -nw $readid $FQDIR/${FQID}_R1.fastq | cut -f1 -d':')
 linenum=$((linenum + 1))
 echo "FASTQ line $linenum"

 sedcommand="$sedcommand -e ${linenum}s/${data}/${datamod}/"
done <<< "$READS"

############# execute the actual fastq editing
echo "sed -i ${sedcommand} ${FQDIR}/${FQID}_R1.fastq"
sed -i ${sedcommand} ${FQDIR}/${FQID}_R1.fastq




sedcommand=""
########### reads came from R2 and have been revcomped in the BAM (i.e. from reverse strand)
READS=$(samtools view -h -f 128 -f 16 $BAMDIR/${BAMID}.dedup.bam Chr01:${SNPPOS}-${SNPPOS} | grep -v "^@" \
 | awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')


echo "$READS"

while read -r snppos offset base readid readstart data throwaway; do
 echo "read to edit: " $snppos $offset $base $readid $f5 $data;
 # search fastq file of that sample for the read ID

 echo "original:" $data
 datarc=$(echo $data | grep '^[ATCG]' - | rev | tr ATCG TAGC)
 datamod=$(echo $data | sed s/./${SNP}/${offset})                               # replace the specific base at offset with the simulated mutation
 datamod2=$(echo $datamod | grep '^[ATCG]' - | rev | tr ATCG TAGC)               # reverse complement the alignment
 echo "proposed:" $datamod2

 linenum=$(grep -m 1 -nw $readid $FQDIR/${FQID}_R2.fastq | cut -f1 -d':')
 linenum=$((linenum + 1))
 echo "FASTQ line $linenum"

 sedcommand="$sedcommand -e ${linenum}s/${datarc}/${datamod2}/"
done <<< "$READS"


######## reads came from R2 and have NOT been revcomped in the BAM (i.e. from forward strand)
READS=$(samtools view -h -f 128 -F 16 $BAMDIR/${BAMID}.dedup.bam Chr01:${SNPPOS}-${SNPPOS} | grep -v "^@" \
 | awk -v pos=$SNPPOS -v clip=0 'BEGIN {OFS=FS="\t"} ; { c=split($6,cigar,"S"); clip=cigar[1]; if(clip ~ /^[0-9]+$/) clip=clip; else clip=0 };  {n=split($10,a,"")} ; {print pos,(pos-$4)+clip+1, a[(pos-$4)+clip+1], $1, $4, $10}')


echo "$READS"
# NOTE: placing quotes around the variable means that the echo output is AS IS.

while read -r snppos offset base readid readstart data throwaway; do
 echo "read to edit: " $snppos $offset $base $readid $f5 $data;
 # search fastq file of that sample for the read ID

 echo "original:" $data
 datamod=$(echo $data | sed s/./${SNP}/${offset}) # replace the specific base at offset with the simulated mutation
 echo "proposed:" $datamod

 linenum=$(grep -m 1 -nw $readid $FQDIR/${FQID}_R2.fastq | cut -f1 -d':')
 linenum=$((linenum + 1))
 echo "FASTQ line $linenum"

 sedcommand="$sedcommand -e ${linenum}s/${data}/${datamod}/"
done <<< "$READS"


############# execute the actual fastq editing
echo "sed -i ${sedcommand} ${FQDIR}/${FQID}_R2.fastq"
sed -i ${sedcommand} ${FQDIR}/${FQID}_R2.fastq

