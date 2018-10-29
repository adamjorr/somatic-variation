# zypy
WIP
Only slightly related to working with heterozygous sequences and not even that much Python

The First Attempt
=================

DiscoSNP++
----------

`$ run_discoSnp++.sh -r \"$(find ../data/ -name *.fastq | sort)\" -m -T`

`$ run_VCF_creator.sh -G ref.fa -p discoRes_k_31_c_4_D_0_P_1_b_0_withlow_coherent.fa -o discosnp.vcf -w`

`$ bash vcf2tree.sh discosnp.vcf`

GATK
----

`$ bash bowtie2aligner.sh ref.fa ../data/ out.bam`

`$ bash gatkcaller.sh ref.fa out.bam`

`$ bash vcf2tree.sh var-calls.vcf`


With Kmer Correction and Trimming
=================================
Run Rcorrector with the script included in the program:
```bash
for F in $(find ../data/ -name '*R1*.fastq'); do perl run_rcorrector.pl -1 $F -2 ${F/R1/R2} -k 32 -t 48 -od ../corrected_data/; done
```

Build a graph and filter on estimated coverage using Khmer. Check [this fork](https://github.com/adamjorr/khmer):
```bash
scripts/load-into-counting.py -ksize 32 -T 48 -M 64e9 khmer_count.graph ../corrected_data/*.fq
for F in $(find ../corrected_data/ -name '*R1*.fq'); do sandbox/slice-paired-reads-by-coverage.py -M 40000 khmer_count.graph $F ${F/R1/R2} $(basename $F fq)_sliced.fq $(basename ${F/R1/R2} fq)_sliced.fq $(basename ${F/R1/} fq)_singletons.fq; done
```

For a parallelized version, use:
```bash
parallel -j ${CORES} 'F={}; G={/.}; sandbox/slice-paired-reads-by-coverage.py -M 40000 khmer_count.graph $F ${F/R1/R2} ${G}_sliced.fq ${G/R1/R2}_sliced.fq ${G/R1/}_singletons.fq ::: $(find ../corrected_data/ -name '*R1*.fq')'
```


Align with BWA, extract the unmapped and poorly mapped (-q argument) reads, then realign them with stampy. **WARNING:** stampy seems to crash with samtools 1.3, but 1.2 seems to work fine.
```bash
bash bwa_aligner.sh ref.fa ../corrected_data/ out1.bam
samtools sort -@ 48 -n -m 2G out1.bam -o out1_sorted.bam
samtools view -@ 48 -b -h -q 13 -f 2 -o out1_mapped.bam -U out1_unmapped.bam out1_sorted.bam
bash stampy_realigner.sh ref.fa out1_unmapped.bam out1_remapped.bam
samtools merge -@ 48 -n -c -p out1_merged.bam out1_mapped.bam out1_remapped.bam
samtools sort -@ 48 -m 2G -o out1_sorted.bam -O bam -T tmp out1_merged.bam
```

Create a new reference. See the instructions [here](https://github.com/samtools/bcftools/wiki/HOWTOs#consensus-calling).
We don't use the vcfutils.pl method because it uses iupac ambiguity codes.
```bash
samtools mpileup -uf ref.fa out1_sorted.bam | bcftools call --threads 48 -mv -Ou | bcftools filter --threads 48 -Oz -o bcftools_calls.vcf.gz -e 'AC==1'
tabix bcftools_calls.vcf.gz
cat ref.fa | bcftools consensus bcftools_calls.vcf.gz -c consensus.chain > consensus.fa
```

Align to the new reference and call variants with GATK.
```bash
bash bwa_aligner.sh consensus.fa ../data/ out2.bam
bash gatkcaller.sh consensus.fa bwa_2.bam
```

This will create a file called var-calls.vcf which will contain the variants. This can be fed into `vcf2tree.sh` to construct trees using various filtering methods.

Deviations
----------

To make things easier to handle in Rcorrector, we concatenated all the reads into two files (one for each pair). This ended up making things a lot more difficult and made it impossible to add read groups to the initial alignment. This *shouldn't* make a difference, but may have affected the consensus reference. I would **not** recommend doing this in the future.


Recent Additions and Simplifications
====================================
Do kmer correction with:
```bash
bash clean_reads.sh ../data/
```

Then align to a reference with bwa:
```bash
bash bwa_aligner.sh ref.fa ./cleaned_reads/sliced/ bwa1.bam
```

Next, realign with Stampy:
```bash
bash stampy_realigner.sh ref.fa bwa1.bam stampy1.bam
```

Now obtain a consensus with


TBC

Average Number of Mismatches per Read
=====================================
AVGMISMATCHES=$(samtools view -h -s .01 ../e_mel_1/e_mel_bwa_1.bam |
samtools view -h -b -u -q ${QUAL} -F4 -o /dev/null -U /dev/stdout |
samtools calmd -e - ../ref.fa |
cut -f10 |
head -n 20 |
tr -d = |
awk '{sum += length} END {printf "%0.0f\n" ,sum / NR}')
