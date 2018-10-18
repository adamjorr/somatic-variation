
all: chr_plot.pdf

genes.txt : ../liftover/e_mel_3_genes.gff3
	cat $< | sed -nE '/^scaffold_[0-9]\t|scaffold_1[01]\t/p' | awk '$$3 == "gene" {print}' | cut -f1,4,5,9 | awk 'BEGIN{OFS="\t";print "chr","start","end","name","gieStain"}{gsub(/.+[nN]ame=/,"",$$4); print $$0, "gvar"}' >$@

chromosomes.txt : ../filtered_bed_excluded.vcf
	bcftools view -h $< | sed -n -E '/ID=scaffold_[0-9],|ID=scaffold_1[01],/p' | tr -d '<>' | sed -nE 's/length=//gp' | cut -d= -f3 | tr , '\t' | awk 'BEGIN{OFS="\t";print "chr","start","end"}{print $$1, 1, $$2}' >$@

positions_with_repeats.txt : ../filtered.vcf
	cat $< | filt_with_replicates.pl -s -g 3 | diploidify.py -v -t vcf | bcftools view -H | cut -f1,2 > $@

positions_no_repeats.txt : ../filtered_bed_excluded.vcf
	cat $< | filt_with_replicates.pl -s -g 3 | diploidify.py -v -t vcf | bcftools view -H | cut -f1,2 > $@

chr_plot.pdf : positions_no_repeats.txt positions_with_repeats.txt chromosomes.txt genes.txt
	Rscript chr_plot.R
