SCRIPTDIR=../scripts
ANALYSISDIR=../analysis
CALLABILITYDIR=$(ANALYSISDIR)/false_negative_rate
DNGDIR=$(ANALYSISDIR)/dng
EMEL3DIR=$(ANALYSISDIR)/e_mel_3

dng_callability.txt : $(CALLABILITYDIR)/dng/results.txt
	cp $< > $@

dng_variants.tsv : $(ANALYSISDIR)/variant_table/var_calls.tsv
	cp $< > $@

filtered_bed_excluded.fa : $(ANALYSISDIR)/gatk4/replicate-filter-only-variable.fa
	cp $< > $@

genomecov.txt : $(ANALYSISDIR)/liftover/e_mel_3_genes.bed $(ANALYSISDIR)/e_mel_3/e_mel_3.bed
	bedtools genomecov -max 1 -i $< -g $(word 2, $^) > $@

MC_variants.bed : dng_variants.tsv
	cat $< | awk 'BEGIN{OFS="\t";} ($53 == "MC"){print $1, $2 - 1, $2}' > $@

MC_closest_features.bed : MC_variants.bed $(ANALYSISDIR)/liftover/e_mel_3_genes.bed
	bedtools closest -a $< -b $(word 2, $^) > $@

MC_closest_genes.txt : MC_closest_features.bed
	awk 'BEGIN{OFS="\t"}; ($11 == "gene"){print $1, $2, $5, $6, $7, $13 }' $< > $@

MC_genes_phytozome_descriptions.tsv : MC_closest_genes.txt
	cat $< | $(SCRIPTDIR)/phytozome_query.py > $@

