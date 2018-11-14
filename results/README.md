```
cat dng_variants.tsv | awk 'BEGIN{OFS="\t";} ($53 == "MC"){print $1, $2 - 1, $2 - 1}' > MC_variants.bed
sort -n -k1,10 -k2,2 e_mel_3_genes.bed > e_mel_3_genes.sorted.bed
bedtools closest -g /storage/2017-03-15-eucalyptus/e_mel_3/e_mel_3.bed -a MC_variants.bed -b /storage/2017-03-15-eucalyptus/liftover/e_mel_3_genes.sorted.bed > MC_closest_features.bed
awk 'BEGIN{OFS="\t"}; ($11 == "gene"){print $1, $2, $5, $6, $7, $13 }' MC_closest_features.bed > MC_closest_genes.txt
awk 'BEGIN{split("scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11",tmp); for (i in tmp) targets[tmp[i]]} $1 in targets {print $0}' e_mel_3_genes.resort.bed > e_mel_3_genes_first11.bed
bedtools genomecov -max 1 -i e_mel_3_genes.resort.bed -g <(head -n11 /storage/2017-03-15-eucalyptus/e_mel_3/e_mel_3.bedtools.genome) > genomecov.txt
```
 * `MC_variants.bed` : BED file of positions of variants on the resistant branch
 * `MC_closest_features.bed` : BED file with annotations closest to the positions mutated on the resistant branch
 * `MC_closest_genes.txt` : Just genes in `MC_closest_featutes.bed` and fewer columns
 * `MC_genes_phytozome_descriptions.tsv` : Manually queried phytozome for genes in `MC_closest_genes.txt`
