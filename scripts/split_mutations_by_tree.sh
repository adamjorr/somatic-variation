#!/bin/bash

#SHUFFLED=$(shuf /short/xf1/Emel/induce_mutations)


# parse command line inputs (bed file and num lines per chunk)
while getopts i:o:n: option
do
 case "${option}"
 in
 i) INFILE=${OPTARG};;
 o) OUTDIR=${OPTARG};;
 n) NLINES=${OPTARG};;
 :) echo "Option -"$OPTARG" requires an argument" >&2;   exit 1  ;;
 *) echo "invalid option $OPTARG" >&2;   exit    ;;
 esac
done

if [ $OPTIND -eq 1 ]; then echo "No options were passed. We really must have options. We simply can't go on without options. Please...options"; exit; fi

shuf $INFILE | split -l $NLINES - ${OUTDIR}/x

cat ${OUTDIR}/xaa >> ${OUTDIR}/mut_M1.txt
cat ${OUTDIR}/xaa >> ${OUTDIR}/mut_M2.txt
cat ${OUTDIR}/xaa >> ${OUTDIR}/mut_M3.txt
cat ${OUTDIR}/xaa >> ${OUTDIR}/mut_M4.txt

cat ${OUTDIR}/xab >> ${OUTDIR}/mut_M1.txt
cat ${OUTDIR}/xab >> ${OUTDIR}/mut_M2.txt
cat ${OUTDIR}/xab >> ${OUTDIR}/mut_M3.txt

cat ${OUTDIR}/xac >> ${OUTDIR}/mut_M1.txt
cat ${OUTDIR}/xac >> ${OUTDIR}/mut_M2.txt

cat ${OUTDIR}/xad >> ${OUTDIR}/mut_M1.txt

cat ${OUTDIR}/xae >> ${OUTDIR}/mut_M2.txt

cat ${OUTDIR}/xaf >> ${OUTDIR}/mut_M3.txt

cat ${OUTDIR}/xag >> ${OUTDIR}/mut_M4.txt

cat ${OUTDIR}/xah >> ${OUTDIR}/mut_M5.txt
cat ${OUTDIR}/xah >> ${OUTDIR}/mut_M6.txt
cat ${OUTDIR}/xah >> ${OUTDIR}/mut_M7.txt
cat ${OUTDIR}/xah >> ${OUTDIR}/mut_M8.txt

cat ${OUTDIR}/xai >> ${OUTDIR}/mut_M6.txt
cat ${OUTDIR}/xai >> ${OUTDIR}/mut_M7.txt
cat ${OUTDIR}/xai >> ${OUTDIR}/mut_M8.txt

cat ${OUTDIR}/xaj >> ${OUTDIR}/mut_M6.txt	# This was incorrect in first run. It had M8,M7!
cat ${OUTDIR}/xaj >> ${OUTDIR}/mut_M7.txt

cat ${OUTDIR}/xak >> ${OUTDIR}/mut_M8.txt

cat ${OUTDIR}/xal >> ${OUTDIR}/mut_M7.txt

cat ${OUTDIR}/xam >> ${OUTDIR}/mut_M6.txt

cat ${OUTDIR}/xan >> ${OUTDIR}/mut_M5.txt



#rm ${TMPDIR}/x*
