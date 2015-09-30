# zypy
WIP
Only slightly related to working with heterozygous sequences and not even that much Python

Example
=======
`$ run_discoSnp++.sh -r \"$(find ../data/ -name *.fastq | sort)\" -m -T`

`$ run_VCF_creator.sh -G ref.fa -p discoRes_k_31_c_4_D_0_P_1_b_0_withlow_coherent.fa -o discosnp.vcf -w`

`$ bash vcf2tree.sh discosnp.vcf`
