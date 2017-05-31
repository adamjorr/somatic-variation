#!/usr/bin/env Rscript
#Rscript plot_tree.R tree.nwk out_file.pdf

library('ape')
library('phytools')
mycolors <- c('lightgreen','yellow','purple','red')

args <- commandArgs(trailingOnly = TRUE)
tree <- read.tree(args[1])
pdf(args[2],8,8)

#Fix labels
tree$tip.label = unlist(lapply(tree$tip.label, FUN=function(x) strsplit(x,'_')[[1]][1]))
m <- regexpr("[0-9]+",tree$tip.label)
tree$tip.label = regmatches(tree$tip.label,m)

plot(tree, no.margin = TRUE, cex = 2)

q()