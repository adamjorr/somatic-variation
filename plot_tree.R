#!usr/bin/env Rscript
#Rscript plot_tree.R tree.nwk out_file.pdf

library(ape)
library(phytools)
mycolors <- c('lightgreen','yellow','purple','red')

args <- commandArgs(trailingOnly = TRUE)
tree <- read.tree(args[1])
pdf(args[2],8,8)
plot(tree, show.tip.label = FALSE, no.margin = TRUE, y.lim = c(0,9))
tree$tip.label = 
tiplabels(tree$tip.label,frame="none")

tree$tip.label = unlist(lapply(tree$tip.label, FUN=function(x) strsplit(x,'_')[[1]][1]))
tree$tip.label = gsub("M2c","M3a",tree$tip.label) #fix an oopsie
tree$tip.label = grep
