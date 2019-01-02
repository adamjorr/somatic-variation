# Variant annotations
```
cat dng_variants.tsv | awk 'BEGIN{OFS="\t";} ($53 == "MC"){print $1, $2 - 1, $2 - 1}' > MC_variants.bed
bedtools sort -i e_mel_3_genes.bed > e_mel_3_genes.sorted.bed 
bedtools closest -a MC_variants.sorted.bed -b e_mel_3_genes_first11.sorted.bed > MC_closest_features.bed
awk 'BEGIN{OFS="\t"}; ($11 == "gene"){print $1, $2, $5, $6, $7, $13 }' MC_closest_features.bed > MC_closest_genes.txt
awk 'BEGIN{split("scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11",tmp); for (i in tmp) targets[tmp[i]]} $1 in targets {print $0}' e_mel_3_genes.resort.bed > e_mel_3_genes_first11.bed
bedtools genomecov -max 1 -i e_mel_3_genes.resort.bed -g <(head -n11 /storage/2017-03-15-eucalyptus/e_mel_3/e_mel_3.bedtools.genome) > genomecov.txt
```
 * `MC_variants.bed` : BED file of positions of variants on the resistant branch
 * `MC_closest_features.bed` : BED file with annotations closest to the positions mutated on the resistant branch
 * `MC_closest_genes.txt` : Just genes in `MC_closest_featutes.bed` and fewer columns
 * `MC_genes_phytozome_descriptions.tsv` : Phytozome queried for genes in `MC_closest_genes.txt`

# Number of unmapped reads
```bash
find . -name 'e_mel_*.bam' -or -wholename './alignment.bam' | parallel --tag samtools view -c -f 4 {}
```

 * ./e_mel_3/e_mel_3.bam   59110579                          
 * ./e_mel_2/e_mel_2.bam   61278970                 
 * ./e_mel_1/e_mel_1.bam   66949528                                                                                                           
 * ./alignment.bam 85187119

# Number of Q=0 reads
```bash
for i in $( find . -name 'e_mel_*.bam' -or -wholename './alignment.bam' );
	do echo $i $(samtools view -F 4 $i | awk '$5=="0" { counter++ } END { print counter }' );
done
```
 * ./alignment.bam 305897775
 * ./e_mel_3/e_mel_3.bam 311941189
 * ./e_mel_2/e_mel_2.bam 349608346
 * ./e_mel_1/e_mel_1.bam 311203429

