library(phangorn)
library(phytools)
library(ggplot2)
library(reshape2)
library(adephylo)



# input files
# must exist prior to running this script
alignment = "~/Documents/github/somatic-variation/positive_control_analyses/analysis0_alignments/aln_0gaps.fa"
true.tree = "~/Documents/github/somatic-variation/positive_control_analyses/tree_em1_physical_structure.phy"
bl.meters = "~/Documents/github/somatic-variation/positive_control_analyses/distances_em1.csv"

# make output directory
outdir = "~/Documents/github/somatic-variation/positive_control_analyses/analysis2_compare_brlens/"
dir.create(outdir)

################ Analysis 3 #####################
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


# merge the physical and the parsimony, because these are on the same tree
physical = data.frame(true.tree.physical$edge, 'physical' = true.tree.physical$edge.length)
parsimony = data.frame(true.tree.parsimony$edge, 'parsimony' = true.tree.parsimony$edge.length)

edges = merge(physical, parsimony, by = c('X1', 'X2'))

p = root(true.tree.parsimony, resolve.root = T, outgroup = "M1")

p2 = ggplot(edges, aes(x = physical, y = parsimony)) + 
    geom_smooth(method='lm', fullrange = T) + 
    geom_point(size = 3) +
    xlab("physical branch length (M)") + 
    ylab("genetic branch length (mutations)") + 
    xlim(c(0, 15))


pdf(file.path(outdir, "brlen_correlation_figure_not_forced_through_origin.pdf"), width=7, height=5)
p2
dev.off()


p3 = ggplot(edges, aes(x = physical, y = parsimony)) + 
    geom_smooth(method='lm', formula=y~x-1, fullrange = T) + 
    geom_point(size = 3) +
    xlab("physical branch length (M)") + 
    ylab("genetic branch length (mutations)") + 
    xlim(c(0, 15))


pdf(file.path(outdir, "brlen_correlation_figure_forced_through_origin.pdf"), width=7, height=5)
p3
dev.off()




# regressions forced through the origin
s = summary(lm(edges$physical ~ edges$parsimony -1))
s1 = summary(lm(edges$physical ~ edges$parsimony))

sink(file.path(outdir, "model_fits.txt"))
print(s)
print(s1)
sink()




# mean root-to-tip branchlength of true tree
true.tree.physical = read.tree(true.tree)
mean(distRoot(true.tree.physical))
