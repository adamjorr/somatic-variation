#!/usr/bin/bash

#Script to replicate the Eucalyptus experiment. Takes the location of the directory containing the raw reads as an argument.

READ_DIR = $1
trap "exit 1" ERR
#clean reads
bash clean_reads.sh $READ_DIR

#make first alignment
mkdir e_mel_1 && cd e_mel_1 #do work in a new directory
bash bwa_aligner.sh -r ../e_grandis_f11_cp.fa -o e_mel_bwa_1.bam -i ../cleaned_reads/sliced/ #first aligned the cleaned reads to the 11 chromosomes of the nuclear grandis reference plus melliodora chloroplast
bash stampy_realigner.sh -r ../e_grandis_f11_cp.fa -i e_mel_bwa_1.bam -o e_mel_realigned_1.bam #clean up the alignment with stampy
bash create_consensus.sh -o e_mel_1.fa -r ../e_grandis_f11_cp.fa -i e_mel_realigned_1.bam #make a new consensus reference from the alignment
cd ..

#make second alignment
mkdir e_mel_2 && cd e_mel_2 #do work in a new directory
bash bwa_aligner.sh -r ../e_mel_1/e_mel_1.fa -o e_mel_bwa_2.bam -i ../cleaned_reads/sliced/ #align the cleaned reads to the first version of our melliodora reference
bash stampy_realigner.sh -r ../e_mel_1/e_mel_1.fa -i e_mel_bwa_2.bam -o e_mel_realigned_2.bam #clean up the alignment with stampy
bash create_consensus.sh -o e_mel_2.fa -r ../e_mel_1/e_mel_1.fa -i e_mel_realigned_2.bam #make a new consensus reference from the alignment
cd ..

#make third alignment
mkdir e_mel_3 && cd e_mel_3 #do work in a new directory
bash bwa_aligner.sh -r ../e_mel_2/e_mel_2.fa -o e_mel_bwa_3.bam -i ../cleaned_reads/sliced/ #align the cleaned up reads to the second version of our melliodora reference
bash stampy_realigner.sh -r ../e_mel_2/e_mel_2.fa -i e_mel_bwa_3.bam -o e_mel_realigned_3.bam #clean up the alignment with stampy
bash create_consensus.sh -o e_mel_3.fa -r ../e_mel_2/e_mel_2.fa -i e_mel_realigned_3.bam #make a new consensus reference from the alignment
cd ..

#call variants
bash bwa_aligner.sh -r e_mel_3/e_mel_3.fa -o alignment.bam -i $READ_DIR #align the raw reads to our new pseudo-reference
bash gatkcaller.sh -o var-calls.vcf -r e_mel_3/e_mel_3.fa -i alignment.bam #use GATK to call variants. we use unfiltered reads so we don't mess up assumptions made by GATK.
#outputs var-calls.vcf, the variant calls

#filter variants and make tree
bash vcf2tree.sh -t 12 -g 3 -i var-calls.vcf -o tree.pdf #the fasta alignment of the SNPs is called cleaned.fasta





