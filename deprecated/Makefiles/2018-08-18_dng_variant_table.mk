var_calls.tsv: ../dng/filter/goodsites.vcf ../e_mel_3/e_mel_3.fa ../liftover/e_mel_3_genes.gff3 ../liftover/e_mel_3_repeatmask.gff
	~/bin/zypy/misc/variant_analysis.R $^

chromosome_plot.pdf: ../dng/filter/goodsites.vcf ../e_mel_3/e_mel_3.fa
	~/bin/zypy/plot_chromosomes.R $^
