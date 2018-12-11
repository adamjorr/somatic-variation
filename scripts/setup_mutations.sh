#!/bin/bash

########################################################################################################
# This script takes a list of genomic positions and their genotypes (from a VCF file) and mutates an allele in each.
# If the genotype is homozygous then it mutates the 'ALT' allele
# If the genotype is heterozygous then it randomly chooses the REF or ALT to mutate
# Each mutation is done under the assumption of a 2:1 Ts:Tv ratio.
#
# Input file format is Tab delimited, (no headers):
# CHROM	POS	VCFREF	VCFALT	REF	ALT
# scaffold_1      5093    G       .       G        G
# scaffold_1      11686   A       .       A        A
# scaffold_1      12008   A       .       A        A
# scaffold_1      15101   C       .       C        C
# scaffold_1      21924   C       .       C        C
#######################################################################################################

# These arrays prove a Ts:Tv probability of 2:1 when sampled randomly
# Really what we have is: if(A) then p(A->C)=0.16, p(A->T)=0.16, p(A->G)=0.68, so this is close to it.

Amut=('G' 'G' 'G' 'G' 'T' 'C')
Tmut=('C' 'C' 'C' 'C' 'A' 'G')
Cmut=('T' 'T' 'T' 'T' 'G' 'A')
Gmut=('A' 'A' 'A' 'A' 'T' 'C')

function mutate() {

        local rando=$(($RANDOM % 6))
        #echo $rando

        #$1 is the allele to be mutated, considering Ts:Tv ration of about 2:1
        case "$1" in
                A) echo ${Amut[$rando]}
                ;;
                C) echo ${Cmut[$rando]}
                ;;
                G) echo ${Gmut[$rando]}
                ;;
                T) echo ${Tmut[$rando]}
                ;;
        esac
}


while read -r CHROM POS REF ALT ROOTREF ROOTALT throwaway;
do
	if [ "$ROOTREF" == "$ROOTALT" ]; then	# homozygous site at the tree root
		new=$(mutate "$ROOTALT")
		echo "$CHROM $POS $ROOTREF $ROOTALT $ROOTREF $new" >> $(pwd)/mutfile.txt
	else					# heterozygous site at the tree root, so randomly pick whether to mutate the ref or alt
		if [ $(($RANDOM % 2)) -eq 1 ]; then
		new=$(mutate "$ROOTREF")
		echo "$CHROM $POS $ROOTREF $ROOTALT $new $ROOTALT" >> $(pwd)/mutfile.txt
		else
		new=$(mutate "$ROOTALT")
		echo "$CHROM $POS $ROOTREF $ROOTALT $ROOTREF $new" >> $(pwd)/mutfile.txt
		fi
	fi
done < /dev/stdin

