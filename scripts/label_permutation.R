#!/usr/bin/env Rscript
# this script produces permutations of samples such that
# the new permutation shares no splits in common with the original
# and replicates are are permuted at random, as long as 
# they satisfy the condition that no one set of 3 replicates can
# have >1 replicate from a single sample

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

ids = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
reps = c("a", "b", "c")
labels = expand.grid(reps, ids)
names(labels) = c("original.rep", "original.branch")
labels$original.id = paste(labels[,2], labels[,1], sep="")

treestr = "(((((M1a:0,M1b:0,M1c:0)M1:698,(M2a:0,M2b:0,M2c:0)M2:638)MD:75,(M3a:0,M3b:0,M3c:0)M3:796)MC:773,(M4a:0,M4b:0,M4c:0)M4:844)MB:372,((((M7a:0,M7b:0,M7c:0)M7:978,(M6a:0,M6b:0,M6c:0)M6:928)MG:111,(M8a:0,M8b:0,M8c:0)M8:1029)MF:112,(M5a:0,M5b:0,M5c:0)M5:1297)ME:358)MA:0;"

tree = read.tree(text=str_replace_all(treestr, "\\([^()]+\\)",""))

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

new.tree.str = write.tree(new.tree)

for(i in ids) {
    labs = new.labels %>% filter(original.branch == i) %>% pull(new.id)
    s = str_c(labs,":0",collapse=",")
    pattern = sprintf("%s:",i)
    replacement = sprintf("(%s)%s:", s, i)
    new.tree.str = str_replace(new.tree.str, pattern, replacement)
}

cat(new.tree.str, "\n")
