library(phangorn)
library(phytools)
library(ggplot2)

alignment = "~/Documents//github/somatic-variation/results/filtered_bed_excluded_no_s.fa"
iqtree = "/Users/roblanfear/Documents/github/somatic-variation/variant_analyses/IQtree/iqtree" # mac version by default

# make alignments with zero to 5 gaps per column, and run IQ tree on each

# load the gapped alignment
EM1 = read.phyDat(alignment, format="fasta", type="DNA")

# Get all alignments with different numbers of gaps, e.g. aln4 allows for 4 gaps per column
# write them to file
# run IQ tree with default settings, plus 1000 bootstraps

aln0 = phyDat(del.colgapsonly(EM1,  threshold = 1/8))
aln0f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_iqtree/aln0.phy"
write.phyDat(aln0, file = aln0f, format = "phylip")
c1 = paste(iqtree, "-s", aln0f, "-b 250", sep=" ")
system(c1)

aln1 = phyDat(del.colgapsonly(EM1,  threshold = 2/8))
aln1f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_iqtree/aln1.phy"
write.phyDat(aln1, file = aln1f, format = "phylip")
c1 = paste(iqtree, "-s", aln1f, "-b 250", sep=" ")
system(c1)

aln2 = phyDat(del.colgapsonly(EM1,  threshold = 3/8))
aln2f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_iqtree/aln2.phy"
write.phyDat(aln2, file = aln2f, format = "phylip")
c2 = paste(iqtree, "-s", aln2f, "-b 250", sep=" ")
system(c2)

aln3 = phyDat(del.colgapsonly(EM1,  threshold = 4/8))
aln3f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_iqtree/aln3.phy"
write.phyDat(aln3, file = aln3f, format = "phylip")
c3 = paste(iqtree, "-s", aln3f, "-b 250", sep=" ")
system(c3)

aln4 = phyDat(del.colgapsonly(EM1,  threshold = 5/8))
aln4f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_iqtree/aln4.phy"
write.phyDat(aln4, file = aln4f, format = "phylip")
c4 = paste(iqtree, "-s", aln4f, "-b 250", sep=" ")
system(c4)

aln5 = phyDat(del.colgapsonly(EM1,  threshold = 6/8))
aln5f = "~/Documents/github/somatic-variation/variant_analyses/analysis0_iqtree/aln5.phy"
write.phyDat(aln5, file = aln5f, format = "phylip")
c5 = paste(iqtree, "-s", aln5f, "-b 250", sep=" ")
system(c5)    