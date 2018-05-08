library(phangorn)
library(phytools)
library(ggplot2)




# input files
# must exist prior to running this script
alignment = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/aln_0gaps.fa"
true.tree = "~/Documents/github/somatic-variation/variant_analyses/tree_em1_physical_structure.phy"
bl.meters = "~/Documents/github/somatic-variation/variant_analyses/distances_em1.csv"
detection = "~/Documents/github/somatic-variation/variant_analyses/analysis2_detection_rate/detection_rate_per_branch.tsv"

# make output directory
outdir = "~/Documents/github/somatic-variation/variant_analyses/analysis3_compare_brlens/"
dir.create(outdir)



################ Analysis 3.1 #####################
# compare branch lengths of mutations vs. meters  #

# load the alignment
EM1 = read.phyDat(alignment, format="fasta", type="DNA")

# load the physical tree structure
true.tree.physical = read.tree(true.tree)

# important to unroot it, because we have no way of polarising the lengths of the internal root 
# branch for the molecular data
true.tree.physical = unroot(true.tree.physical)
# if you want parsimony branch lengths
true.tree.genomic = acctran(true.tree.physical, EM1)

physical = data.frame(true.tree.physical$edge, 'physical' = true.tree.physical$edge.length)
genomic = data.frame(true.tree.genomic$edge, 'genomic' = true.tree.genomic$edge.length)
edges = merge(physical, genomic, by = c('X1', 'X2'))

# now we label the branches by reading in the dataframe of labelled branch lengths
bl = read.csv(bl.meters)
bl$physical = bl$distance_cm/100
bl = bl[,c('branch', 'physical')]
bl$branch = as.character(bl$branch)

# merge the two root branches on bl
new = c("B->E", bl$physical[which(bl$branch=="A->E")] + bl$physical[which(bl$branch=="A->B")])
bl = rbind(bl, new)

# this is to facilitate the merge, which otherwise appears to hang
bl$physical = as.character(as.numeric(bl$physical))
edges$physical = as.character(edges$physical)

# merge the data frame
edges = merge(edges, bl, by.y='physical')
edges$physical = as.numeric(edges$physical)

edges = edges[,-c(2,3)]

# now we add the correction factors from analysis 2
dr = read.delim(detection, stringsAsFactors = F)

# create the branch that merges across the root
newdr = c("B->E", 
          dr$recovered[which(dr$branch=="A->E")] + dr$recovered[which(dr$branch=="A->B")],
          dr$not.recovered[which(dr$branch=="A->E")] + dr$not.recovered[which(dr$branch=="A->B")],
          NA,
          'internal')
        
dr = rbind(dr, newdr)
dr$recovered = as.numeric(dr$recovered)
dr$not.recovered = as.numeric(dr$not.recovered)
dr$recovery.rate = dr$recovered / (dr$recovered + dr$not.recovered)
dr = dr[-c(2,3)]

edges = merge(edges, dr, by = 'branch')

edges$genomic.corrected = edges$genomic / edges$recovery.rate

ggplot(edges, aes(x = physical, y = genomic)) + 
    geom_smooth(method='lm') + 
    geom_point(aes(colour = branch.type), size = 3) + 
    xlab("physical branch length (M)") + 
    ylab("genetic branch length (mutations)")

ggplot(edges, aes(x = physical, y = genomic.corrected)) + 
    geom_smooth(method='lm') + 
    geom_point(aes(colour = branch.type), size = 3) + 
    xlab("physical branch length (M)") + 
    ylab("genetic branch length (mutations)")


# regressions, forced through the origin and then not
summary(lm(edges$physical ~ edges$genomic -1))
summary(lm(edges$physical ~ edges$genomic))

summary(lm(edges$physical ~ edges$genomic.corrected -1))
summary(lm(edges$physical ~ edges$genomic.corrected))

################ Analysis 3.2 #####################
# Mutation rate calculations                      #

# average mutations per meter
total.mutations = sum(edges$genomic)
total.meters = sum(edges$physical)

detection.rate = read.delim(detection, stringsAsFactors = F)
mean.detection = sum(detection.rate$recovered) / (sum(detection.rate$recovered) + sum(detection.rate$not.recovered))

estimated.mutations = total.mutations / mean.detection

mutations.per.meter = estimated.mutations / total.meters
