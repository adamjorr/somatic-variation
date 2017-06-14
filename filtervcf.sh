#!/usr/bin/bash

input=$1
TMPDIR=$(mktemp -d vcf2tree_tmp_XXXXXXXX)
FILTERTAB=$(mktemp --tmpdir=$TMPDIR --suffix=.tab tmp_filtered_XXXXXXXX)
SITESFA=$(mktemp --tmpdir=$TMPDIR --suffix=.fasta tmp_sites_XXXXXXXX)
GTAB=$(mktemp --tmpdir=$TMPDIR --suffix=.tab tmp_gt_XXXXXXXX)
trap "exit 1" ERR
trap 'rm -rf $TMPDIR' EXIT INT TERM HUP


bcftools filter -g 50 -i 'DP <= 500 && ExcessHet <= 40' $input |
bcftools view -m2 -M2 -v snps |
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AF[\t%TGT]\n' > $FILTERTAB

bcftools query -l $input | tr '\n' '\t' > $GTAB
echo >> $GTAB
cat $FILTERTAB | cut -f 7- >> $GTAB
cat $GTAB | tr -d '/' | datamash transpose | awk '{print ">"$1; for (i=2; i<NF; i++) printf $i; print "\n"'} > $SITESFA

# raxmlHPC-PTHREADS-AVX2 -s $SITESFA -T 4 -m GTRCAT -p 1234 -n sites
raxmlHPC -T 4 -f a -s $SITESFA -n nwk -m ASC_GTRGAMMA -w $TMPDIR --asc-corr=lewis -p 12345 -x 12345 -# 100
cat $TMPDIR/RAxML_bestTree.nwk >$2
