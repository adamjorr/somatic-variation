library(phangorn)
library(phytools)

outdir = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/"
dir.create(outdir)

# input files
# must exist prior to running this script
original = "~/Documents//github/somatic-variation/results/filtered_bed_excluded_no_s.fa"

################ Analysis 0 ########################
# create alignments with different numbers of gaps #

# copy the alignment 
system(paste("cp", original, "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/original.fa"))
# make n's into gaps
system("sed -i '' 's/N/-/g' ~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/original.fa")
# remove lowercase 'a' to change e.g. M1a -> M1
system("sed -i '' 's/a/''/g' ~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/original.fa")

# load the alignment with gaps
alignment = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/original.fa"


# make alignments with zero to 5 gaps per column, from the original which has up to 7

# load the gapped alignment
EM1 = read.phyDat(alignment, format="fasta", type="DNA")

# Get all alignments with different numbers of gaps, e.g. aln4 allows for 4 gaps per column
# write them to file
aln0 = phyDat(del.colgapsonly(EM1,  threshold = 1/8))
aln0f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/aln_0gaps.fa"
write.phyDat(aln0, file = aln0f, format = "fasta")

aln1 = phyDat(del.colgapsonly(EM1,  threshold = 2/8))
aln1f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/aln_1gap.fa"
write.phyDat(aln1, file = aln1f, format = "fasta")

aln2 = phyDat(del.colgapsonly(EM1,  threshold = 3/8))
aln2f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/aln_2gaps.fa"
write.phyDat(aln2, file = aln2f, format = "fasta")

aln3 = phyDat(del.colgapsonly(EM1,  threshold = 4/8))
aln3f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/aln_3gaps.fa"
write.phyDat(aln3, file = aln3f, format = "fasta")

aln4 = phyDat(del.colgapsonly(EM1,  threshold = 5/8))
aln4f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/aln_4gaps.fa"
write.phyDat(aln4, file = aln4f, format = "fasta")

aln5 = phyDat(del.colgapsonly(EM1,  threshold = 6/8))
aln5f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/aln_5gaps.fa"
write.phyDat(aln5, file = aln5f, format = "fasta")
