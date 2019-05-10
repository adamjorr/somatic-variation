filtered_including_repeats.vcf: ../filtered.vcf
	cat $< | filt_with_replicates.pl -s -g 3 | diploidify.py -v -t vcf > $@

var_calls.tsv: filtered_including_repeats.vcf ../e_mel_3/e_mel_3.fa ../liftover/e_mel_3_genes.gff3 ../liftover/e_mel_3_repeatmask.gff
	~/bin/zypy/misc/variant_analysis.R $^
