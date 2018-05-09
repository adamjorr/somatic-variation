library(phangorn)
library(phytools)
library(ggplot2)
library(reshape2)




# input files
# must exist prior to running this script
alignment = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/aln_0gaps.fa"
true.tree = "~/Documents/github/somatic-variation/variant_analyses/tree_em1_physical_structure.phy"
bl.meters = "~/Documents/github/somatic-variation/variant_analyses/distances_em1.csv"
detection = "~/Documents/github/somatic-variation/variant_analyses/analysis2_detection_rate/detection_rate_per_branch.tsv"
iqtree = "~/Documents/github/somatic-variation/variant_analyses/iqtree" # mac version by default


# make output directory
outdir = "~/Documents/github/somatic-variation/variant_analyses/analysis3_compare_brlens/"
dir.create(outdir)

# estimated diploid genome size for E. melliodora in bp
genome.size = 1000000000 


################ Analysis 3.1 #####################
# compare branch lengths of mutations vs. meters  #

# load the alignment
EM1 = read.phyDat(alignment, format="fasta", type="DNA")

# load the physical tree structure
true.tree.physical = read.tree(true.tree)

# important to unroot it, because we have no way of polarising the lengths of the internal root 
# branch for the molecular data
true.tree.physical = unroot(true.tree.physical)

# get parsimony branch lengths
true.tree.parsimony = acctran(true.tree.physical, EM1)

# run IQtree assuming the true tree, to get ML branch lenths
system(paste("cp", alignment, file.path(outdir, "aln_0gaps.fa"))) # copy alignment
system(paste("cp", true.tree, outdir)) # copy tree
command = paste(iqtree, "-s", file.path(outdir, "aln_0gaps.fa"), "-redo -m K2P -g ", true.tree, sep=" ")
system(command)
ml.tree.file = file.path(outdir, "aln_0gaps.fa.treefile")
true.tree.ml = read.tree(ml.tree.file)
# update edge lengths in ML tree to be substitutions, not substitutions per site
true.tree.ml$edge.length = true.tree.ml$edge.length * length(as.character(EM1)[1,])


# merge the physical and the parsimony, because these are on the same tree
physical = data.frame(true.tree.physical$edge, 'physical' = true.tree.physical$edge.length)
parsimony = data.frame(true.tree.parsimony$edge, 'parsimony' = true.tree.parsimony$edge.length)
ml = data.frame(true.tree.ml$edge, 'ml' = true.tree.ml$edge.length)

edges = merge(physical, parsimony, by = c('X1', 'X2'))
edges$ml = NA

p = root(true.tree.parsimony, resolve.root = T, outgroup = "M1")
m = root(true.tree.ml, resolve.root = T, outgroup = "M1")
# a long-winded way to match up branch lengths
for(i in 1:14){
    # i is the number in the X2 column of edges
    # for each i, we need to figure out the edge in the ML tree
    if(i %in% edges$X2){ # one branch (the root) isn't
        p.tips = Descendants(p, i, type = "tips")
        p.tips = p$tip.label[p.tips[[1]]]
        print("new")
        print(p.tips)
        # now find the corresponding ML branch for this
        for(k in 1:14){
            
            m.tips = Descendants(m, k, type = "tips")
            m.tips = m$tip.label[m.tips[[1]]]
            if(length(setdiff(p.tips, m.tips)) + length(setdiff(m.tips, p.tips)) == 0){
                print(m.tips)
                
                # now k is the number in the X2 column of ml matrix
                ml.edge.length = ml$ml[which(ml$X2 == k)]
                edges$ml[which(edges$X2==i)] = ml.edge.length
                break
            }
        }
        
    }
}
    


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

edges$parsimony.corrected = edges$parsimony / edges$recovery.rate
edges$ml.corrected = edges$ml / edges$recovery.rate

edges.melt = melt(edges, id.vars = c("branch", "physical", "recovery.rate", "branch.type"))


ggplot(edges.melt, aes(x = physical, y = value)) + 
    geom_smooth(method='lm') + 
    geom_point(aes(colour = branch.type), size = 3) +
    facet_wrap(~variable, scales = "free_y") +
    xlab("physical branch length (M)") + 
    ylab("genetic branch length (mutations)")


# regressions, forced through the origin and then not
summary(lm(edges$physical ~ edges$parsimony -1))
summary(lm(edges$physical ~ edges$parsimony))
summary(lm(edges$physical ~ edges$parsimony.corrected -1))
summary(lm(edges$physical ~ edges$parsimony.corrected))
summary(lm(edges$physical ~ edges$ml -1))
summary(lm(edges$physical ~ edges$ml))
summary(lm(edges$physical ~ edges$ml.corrected -1))
summary(lm(edges$physical ~ edges$ml.corrected))

################ Analysis 3.2 #####################
# Mutation rate calculations                      #

# average mutations per meter
total.mutations.ml = sum(edges$ml)
total.mutations.p = sum(edges$parsimony)

total.meters = sum(edges$physical)

detection.rate = read.delim(detection, stringsAsFactors = F)
mean.detection = sum(detection.rate$recovered) / (sum(detection.rate$recovered) + sum(detection.rate$not.recovered))

estimated.mutations.ml = total.mutations.ml / mean.detection
estimated.mutations.p = total.mutations.p / mean.detection

mutations.per.meter.ml = estimated.mutations.ml / total.meters
mutations.per.meter.p = estimated.mutations.p / total.meters

mutations.per.meter.ml
mutations.per.meter.p

# average mutations per generation
# assuming all mutations are hertiable
# and that a generation is characterised by the root to tip branch length in M

t = read.tree(true.tree)
mean.root.to.tip.meters = mean(diag(vcv.phylo(t)))

mutations.per.generation.ml = mutations.per.meter.ml * mean.root.to.tip.meters
mutations.per.generation.p = mutations.per.meter.p * mean.root.to.tip.meters

mutations.per.generation.ml
mutations.per.generation.p

# average mutations per site per generation

mutations.per.generation.ml / genome.size
mutations.per.generation.p / genome.size
