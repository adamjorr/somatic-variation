library(phangorn)
library(phytools)
library(ggplot2)

outdir = "~/Documents/github/somatic-variation/positive_control_analyses/analysis1_positive_control/"
dir.create(outdir)

# input files
# these files should exist prior to the anlaysis
alignment = "~/Documents/github/somatic-variation/positive_control_analyses/analysis0_alignments/aln_0gaps.fa"
true.tree = "~/Documents/github/somatic-variation/positive_control_analyses/tree_em1_physical_structure.phy"
iqtree = "~/Documents/github/somatic-variation/positive_control_analyses/iqtree" # mac version by default


################ Analysis 1 #####################
# compare the true tree to the inferred tree(s) #

# copy the alignment into this folder, to keep clear what we're doing
system(paste("cp", alignment, file.path(outdir, "aln_0gaps.fa")))

# load the alignment
EM1 = read.phyDat(file.path(outdir, "aln_0gaps.fa"), format="fasta", type="DNA")

# get all possible trees for 8 tips
trees = allTrees(8, tip.label = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"))

# get all possible parsimony scores
scores = lapply(trees, parsimony, EM1)
scores = unlist(scores)

# the best trees have the lowest scores
best = which(scores == min(scores))
best.trees = trees[best]

# All path distances from the true tree
true.tree = read.tree(true.tree)
all.pds = path.dist(trees, true.tree)
all.pds = sort(all.pds)
true.pds = path.dist(best.trees, true.tree)

# p.values: the position of the true pd in the ranked list of observed pds
get_p = function(pd, all.pds){
    return(max(which(all.pds==pd))/length(all.pds))
}

pvals = lapply(true.pds, get_p, all.pds = all.pds)

# Now let's get the maximum likelihood tree for that alignment.
# run IQ-TREE with default settings and 1000 bootstraps
command = paste(iqtree, "-s", file.path(outdir, "aln_0gaps.fa"), "-b 1000 -safe -redo -m JC", sep=" ")
system(command)
ml.tree.file = file.path(outdir, "aln_0gaps.fa.treefile")

ml.tree = read.tree(ml.tree.file)

ml.pd = path.dist(ml.tree, true.tree)
ml.pval = get_p(ml.pd, all.pds)

# what PD would be significant at the 5% level? 
five_percent = all.pds[as.integer(length(all.pds)*0.05)]

# plot the histogram of Path Distances, with lines for the observed distances
# The histogram shows the path distance of all possible 10395 trees of 8 tips from the true tree. T
# The true phylogeny has a path distance of zero. 
# The trees that best explain the alignment according to maximum parsimony are shown as dashed red lines. 
# The solid red line shows the minimum path distance for a tree to be considered significantly closer to the true tree than would be expected by chance.
# These trees are a lot closer to the true tree (nearer zero) than we would expect by chance (solid red line).
pd.df = data.frame(Path.Distance = all.pds)
p1 = ggplot(pd.df, aes(x = Path.Distance)) + 
        geom_histogram(bins = 50) + 
        geom_vline(xintercept = ml.pd, colour = 'blue', size = 1, alpha = 0.5) +
        geom_vline(xintercept = true.pds, colour = 'red', linetype = 'dashed', size = 0.5) +
        geom_vline(xintercept = five_percent, colour = 'red', size = 1, alpha = 0.5) +
        xlab("Path distance of inferred tree from true tree") +
        labs(title = "Positive control analysis")

# now let's output the data and the plot
pdf(file.path(outdir, "positive_control_figure.pdf"), width=7, height=5)
p1
dev.off()

# Write a table of results to best_trees.tsv
r = data.frame(path.distance = true.pds, parsimony.score = scores[best], p.value = unlist(pvals), tree = write.tree(best.trees), type = "maximum parsimony")

r$tree = as.character(r$tree)
r$type = as.character(r$type)

ml = c(ml.pd, NA, ml.pval, write.tree(ml.tree), "maximum likelihood")

r = rbind(r, ml)

write.table(r, file.path(outdir, "best_trees.tsv"), quote=FALSE, sep='\t', col.names = NA)


