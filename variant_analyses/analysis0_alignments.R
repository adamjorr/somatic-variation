library(phangorn)
library(phytools)

outdir = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/"
dir.create(outdir)

# input files
# must exist prior to running this script
original = "~/Documents//github/somatic-variation/results/filtered_bed_excluded.fa"

################ Analysis 0 ########################
# create alignments with no  gaps #

# copy the alignment 
system(paste("cp", original, "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/original.fa"))
# make n's into gaps
system("sed -i '' 's/N/-/g' ~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/original.fa")
# remove lowercase 'a' to change e.g. M1a -> M1
system("sed -i '' 's/a/''/g' ~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/original.fa")

# load the alignment with gaps
alignment = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/original.fa"


# load the gapped alignment
EM1 = read.phyDat(alignment, format="fasta", type="DNA")

# Write alignment with no gaps
# write them to file
aln0 = phyDat(del.colgapsonly(EM1,  threshold = 1/8))
aln0f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/aln_0gaps.fa"
write.phyDat(aln0, file = aln0f, format = "fasta")

