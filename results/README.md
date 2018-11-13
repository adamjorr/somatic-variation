cat dng_variants.tsv | awk 'BEGIN{OFS="\t";} ($53 == "MC"){print $1, $2 - 1, $2 - 1}' > MC_variants.bed
sort -n -k1,10 -k2,2 e_mel_3_genes.bed > e_mel_3_genes.sorted.bed
bedtools closest -g /storage/2017-03-15-eucalyptus/e_mel_3/e_mel_3.bed -a MC_variants.bed -b /storage/2017-03-15-eucalyptus/liftover/e_mel_3_genes.sorted.bed > MC_closest_features.bed
awk 'BEGIN{OFS="\t"}; ($11 == "gene"){print $1, $2, $5, $6, $7, $13 }' MC_closest_features.bed > MC_closest_genes.txt
