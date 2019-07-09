library(phangorn)
library(phytools)

outdir = "~/Documents/github/somatic-variation/positive_control_analyses/analysis0_alignments/"
dir.create(outdir)

# input files
# must exist prior to running this script
original = "~/Documents//github/somatic-variation/results/filtered_bed_excluded.fa"

################ Analysis 0 ########################
# create alignments with no  gaps #

# copy the alignment 
system(paste("cp", original, file.path(outdir, "original.fa")))

# make n's into gaps
system(paste("sed -i '' 's/N/-/g'", file.path(outdir, "original.fa")))

# remove lowercase 'a' to change e.g. M1a -> M1
system(paste("sed -i '' 's/a/''/g'", file.path(outdir, "original.fa")))

# load the alignment with gaps
alignment = file.path(outdir, "original.fa")
EM1 = read.phyDat(alignment, format="fasta", type="DNA")

# Write alignment with no gaps
aln0 = phyDat(del.colgapsonly(EM1,  threshold = 1/8))
aln0f = file.path(outdir, "aln_0gaps.fa")
write.phyDat(aln0, file = aln0f, format = "fasta")
