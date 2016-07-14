#!usr/bin/bash
#bash vcf2tree.sh file.vcf
#Takes file.vcf, filters it using replicate info at various stringencies, and plots the trees.

WHEREAMI=`dirname $0`
PATHTOVCFALNER=$WHEREAMI/vcf_tab_to_fasta_alignment.pl
PATHTOZYPY=$WHEREAMI

#Now we make another directory and do the same thing, but don't allow missing data AND only use variable sites.
mkdir -p strict || exit 1
cd strict || exit 1

perl $PATHTOZYPY/filt_with_replicates.pl -s -g 3 <../$1 >filtered.vcf || exit 1
vcf-to-tab <filtered.vcf >filtered.tab || exit 1
perl $PATHTOVCFALNER -i filtered.tab > cleaned.fasta || exit 1
rm filtered.tab_clean || exit 1
python2 ${PATHTOZYPY}/diploidify.py -i cleaned.fasta -t fasta -o cleaned.dip.phylip-relaxed -p phylip-relaxed -v || exit 1

mkdir -p tree || exit 1
cd tree || exit 1

raxmlHPC -T 4 -f a -s ../cleaned.dip.phylip-relaxed -n nwk -m ASC_GTRGAMMA --asc-corr=lewis -p 12345 -x 12345 -# 100 || exit 1

cd ..

Rscript ${PATHTOZYPY}/plot_tree.R tree/RAxML_bestTree.nwk tree.pdf || exit 1

cd ..


exit
