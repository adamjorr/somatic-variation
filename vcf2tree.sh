#!usr/bin/bash
#bash vcf2tree.sh file.vcf
#Takes file.vcf, filters it using replicate info at various stringencies, and plots the trees.

WHEREAMI=`dirname $0`
PATHTOVCFALNER=$WHEREAMI/vcf_tab_to_fasta_alignment.pl
PATHTOZYPY=$WHEREAMI

#Standard filtering

perl $PATHTOZYPY/filt_with_replicates.pl -g 3 <$1 >filtered.vcf || exit 1
vcf-to-tab <filtered.vcf >filtered.tab || exit 1
perl $PATHTOVCFALNER -i filtered.tab > cleaned.fasta || exit 1
rm filtered.tab_clean || exit 1
python ${PATHTOZYPY}/diploidify.py -i cleaned.fasta -t fasta -o cleaned.dip.phylip-relaxed -p phylip-relaxed || exit 1

#This keeps the ambiguous notation to collect some stats used by supported_sites.py
python ${PATHTOZYPY}/diploidify.py -i cleaned.fasta -t fasta -o cleaned.het.phylip-relaxed -p phylip-relaxed -s || exit 1

#Make the tree in another folder
mkdir -p tree || exit 1
cd tree || exit 1

raxmlHPC -T 4 -f a -s ../cleaned.dip.phylip-relaxed -n nwk -m GTRGAMMA -p 12345 -x 12345 -# 100 || exit 1

cd ..

#Get stats to paint on the tree & do it.
python ${PATHTOZYPY}/supported_sites.py -t tree/RAxML_bipartitions.nwk -a cleaned.het.phylip-relaxed -s -o supp_table.txt -z zygosity.txt -n || exit 1
Rscript ${PATHTOZYPY}/root_homozygote.R zygosity.txt supp_table.txt tree/RAxML_bipartitions.nwk || exit 1

#Now we make another directory and do the same thing, but only use variable sites.
mkdir -p onlyvariable || exit 1
cd onlyvariable || exit 1
python ${PATHTOZYPY}/diploidify.py -i ../cleaned.fasta -t fasta -o cleaned.dip.phylip-relaxed -p phylip-relaxed -v || exit 1
python ${PATHTOZYPY}/diploidify.py -i ../cleaned.fasta -t fasta -o cleaned.het.phylip-relaxed -p phylip-relaxed -v -s || exit 1

mkdir -p tree || exit 1
cd tree || exit 1

raxmlHPC -T 4 -f a -s ../cleaned.dip.phylip-relaxed -n nwk -m ASC_GTRGAMMA --asc-corr=lewis -p 12345 -x 12345 -# 100 || exit 1

cd ..

python ${PATHTOZYPY}/supported_sites.py -t tree/RAxML_bipartitions.nwk -a cleaned.het.phylip-relaxed -s -o supp_table.txt -z zygosity.txt -n || exit 1
Rscript ${PATHTOZYPY}/root_homozygote.R zygosity.txt supp_table.txt tree/RAxML_bipartitions.nwk || exit 1

cd ..

#Now we make another directory and do the same thing, but don't allow missing data AND only use variable sites.
mkdir -p strict || exit 1
cd strict || exit 1

perl $PATHTOZYPY/filt_with_replicates.pl -s -g 3 <../$1 >filtered.vcf || exit 1
vcf-to-tab <filtered.vcf >filtered.tab || exit 1
perl $PATHTOVCFALNER -i filtered.tab > cleaned.fasta || exit 1
rm filtered.tab_clean || exit 1
python ${PATHTOZYPY}/diploidify.py -i cleaned.fasta -t fasta -o cleaned.dip.phylip-relaxed -p phylip-relaxed -v || exit 1
python ${PATHTOZYPY}/diploidify.py -i cleaned.fasta -t fasta -o cleaned.het.phylip-relaxed -p phylip-relaxed -v -s || exit 1

mkdir -p tree || exit 1
cd tree || exit 1

raxmlHPC -T 4 -f a -s ../cleaned.dip.phylip-relaxed -n nwk -m ASC_GTRGAMMA --asc-corr=lewis -p 12345 -x 12345 -# 100 || exit 1

cd ..

python ${PATHTOZYPY}/supported_sites.py -t tree/RAxML_bipartitions.nwk -a cleaned.het.phylip-relaxed -s -o supp_table.txt -z zygosity.txt -n || exit 1
Rscript ${PATHTOZYPY}/root_homozygote.R zygosity.txt supp_table.txt tree/RAxML_bipartitions.nwk || exit 1

cd ..


exit
