library(phangorn)
library(phytools)
library(ggplot2)

alignment = "~/github/somatic-variation/results/filtered_bed_excluded.fa"
ml.treefile = "/Users/roblanfear/Dropbox/Projects_Current/eucalypt_somatic/IQtree/filtered_bed_included.fa.treefile"

################ Analysis 1 #####################
# compare the true tree to the inferred tree(s) #

# load the alignment
EM1 = read.phyDat(alignment, format="fasta", type="DNA")

# get all possible trees for 8 tips
trees = allTrees(8, tip.label = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"))

# here's the parsimony score for a random tree
parsimony(tree = trees[[100]], EM1)

# get all possible parsimony scores
scores = lapply(trees, parsimony, EM1)
scores = unlist(scores)

hist(scores)
min(scores)

best = which(scores == min(scores))
best.trees = trees[best]

# branch lenghts
t1 = acctran(tree = best.trees[[1]], EM1)
plot(t1)
add.scale.bar(length = 100)
t2 = acctran(tree = best.trees[[2]], EM1)

# All path distances from the true tree
true.tree = read.tree('~/Dropbox/trees_sampled/EM1/tree_em1.phy')
all.pds = path.dist(trees, true.tree)
all.pds = sort(all.pds)
true.pds = path.dist(best.trees, true.tree)

# p.values: the position of the true pd in the ranked list of observed pds
pval.1 = max(which(all.pds==true.pds[[1]]))/length(all.pds)
pval.2 = max(which(all.pds==true.pds[[2]]))/length(all.pds)

# what PD would be significant at the 5% level? 
five_percent = all.pds[as.integer(length(all.pds)*0.05)]

# plot the histogram of Path Distances, with lines for the observed distances
pd.df = data.frame(Path.Distance = all.pds)
ggplot(pd.df, aes(x = Path.Distance)) + 
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = true.pds, colour = 'red', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = five_percent, colour = 'red', size = 2, alpha = 0.5)

# build the above plot sequenctially for talks
theme_set(theme_gray(base_size = 20))
ggplot(pd.df, aes(x = Path.Distance)) + 
    geom_histogram(bins = 30) + 
    xlab("Path Distance to structure of physical tree")


ggplot(pd.df, aes(x = Path.Distance)) + 
  geom_histogram(bins = 30) + 
  geom_vline(xintercept = five_percent, colour = 'red', size = 4, alpha = 0.5) +
    xlab("Path Distance to structure of physical tree")


ggplot(pd.df, aes(x = Path.Distance)) + 
  geom_histogram(bins = 30) + 
  geom_vline(xintercept = true.pds, colour = 'red', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = five_percent, colour = 'red', size = 4, alpha = 0.5) +
    xlab("Path Distance to structure of physical tree")



  
# conclusion: the inferred tree is significantly closer to the observed physical topology than can be explained by chance

# we can do the same with the ML tree from IQtree
ml.tree = read.tree(ml.treefile)
ml.tree = root(ml.tree, outgroup = c('M5', 'M6', 'M7', 'M8'), resolve.root = TRUE)
plot(ml.tree)
obj<-cophylo(ml.tree, true.tree)
plot(obj, lwd=3)

association <- cbind(true.tree$tip.label, true.tree$tip.label)
cophyloplot(true.tree, ml.tree, assoc=association, length.line=4, space=28, gap=3)

ml.pds = path.dist(ml.tree, true.tree)
pval.ml = max(which(all.pds==ml.pds))/length(all.pds)

pd.df = data.frame(Path.Distance = all.pds)
ggplot(pd.df, aes(x = Path.Distance)) + 
  geom_histogram(bins = 50) + 
  geom_vline(xintercept = ml.pds, colour = 'red', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = five_percent, colour = 'red', size = 2, alpha = 0.5)


