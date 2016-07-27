#!/bin/sh

input=filtered.vcf

cat $input |
	bcftools filter -g 50 -i 'DP <= 500 && ExcessHet <= 40' |
	bcftools view -m2 -M2 -v snps |
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AF[\t%TGT]\n' > filtered.tab

cat $input | bcftools query -l | tr '\n' '\t' > gt.tab 
echo >> gt.tab
cat filtered.tab | cut -f 7- >> gt.tab
cat gt.tab | tr -d '/' | datamash transpose | awk '{print ">"$1; for (i=2; i<NF; i++) printf $i; print "\n"'} > sites.fasta


rm RAxML_*.sites
raxmlHPC-PTHREADS-AVX2 -s sites.fasta -T 4 -m GTRCAT -p 1234 -n sites
