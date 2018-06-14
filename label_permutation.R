# this script produces permutations of samples such that
# the new permutation shares no splits in common with the original
# and replicates are are permuted at random, as long as 
# they satisfy the condition that no one set of 3 replicates can
# have >1 replicate from a single sample

library(ape)
library(phangorn)
library(dplyr)

ids = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
reps = c("a", "b", "c")
labels = expand.grid(reps, ids)
names(labels) = c("original.rep", "original.branch")
labels$original.id = paste(labels[,2], labels[,1], sep="")

tree = read.tree("~/Documents/github/somatic-variation/variant_analyses/tree_em1_physical_structure.phy")

shuffle.tree = function(original.tree){
    # shuffle tip labels on a tree until the new tree shares NO splits with the old one
    # this is just a tree with the maximum RF dist from the original tree
    
    max.rf = 2*(length(original.tree$tip.label) - 3) # max RF dist
    rf = 0
    
    while(rf < max.rf){
        new.tree = original.tree
        new.tree$tip.label = sample(new.tree$tip.label)
        rf = RF.dist(new.tree, original.tree)
    }
    return(new.tree)
}

shuffle.labels = function(original.labels){
    # shuffle labels so that no more than one branch ID
    # occurs in each new triplet
    # there are better ways than a while loop
    # but it works fine
    
    max.shared = 2
    
    while(max.shared>1){
        labels$new.branch = sample(labels$original.branch)
        labels = labels %>% 
                group_by(original.branch, new.branch) %>% 
                mutate(occ = n())
        max.shared = max(labels$occ)
    }
    
    labels = as.data.frame(labels)
    labels = labels[order(labels$new.branch),]
    labels$new.rep = c("a", "b", "c")
    labels$new.id = paste(labels$new.branch, labels$new.rep, sep = "")
    return(labels)
}

new.labels = shuffle.labels(labels)
new.tree = shuffle.tree(tree)