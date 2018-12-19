SCRIPTDIR=../../scripts/
VARANALYSIS=$(SCRIPTDIR)/variant_analysis.R
CHROMPLOT=$(SCRIPTDIR)/plot_chromosomes.R

var_calls.tsv: ../dng/filter/goodsites.vcf ../e_mel_3/e_mel_3.fa ../e_mel_3/e_mel_3_genes.gff3 ../e_mel_3/e_mel_3_repeatmask.gff
	$(VARANALYSIS) $^

chromosome_plot.pdf: ../dng/filter/goodsites.vcf ../e_mel_3/e_mel_3.fa
	$(CHROMPLOT) $^
