#!/usr/bin/Rscript
#Rscript root_homozygote.R homozygosity_table.txt support_table.txt tree.nwk

#Root the tree with the most homozygous taxon and display some potentially-useful support information.

library(ape)
library(phytools)
mycolors <- c('lightgreen','yellow','purple','red')
pdf('tree.pdf',8,8)

args <- commandArgs(trailingOnly = TRUE)
counts <- read.table(args[1],header = TRUE, sep="\t")
largest <- counts[counts$percent_homozygous==max(counts$percent_homozygous),]
largesample <- largest[[1]]

supportsites <- read.table(args[2],header = TRUE, sep="\t")
splitsupports <- sapply(supportsites[,1],function(x) strsplit(toString(x),','))[[1]]

tree <- read.tree(args[3])
largenodenum <- which(tree$tip.label==largesample)
tree <- reroot(tree,largenodenum,tree$edge.length[which(tree$edge[,2]==largenodenum)])
othersamples <- tree$tip.label[tree$tip.label != largesample]
plot(tree, show.tip.label = FALSE, no.margin = TRUE, y.lim = c(0,9))

numtips <- length(tree$tip.label)
#intedges <- which(!apply(tree$edge,1,function(x) x[1] <= numtips | x[2] <= numtips),arr.ind = TRUE)

for (i in 1:length(supportsites[,1])){
  intaxa <- strsplit(toString(supportsites[i,1]),',')[[1]]
  outtaxa <- setdiff(tree$tip.label,intaxa)
  
  #intaxaedges <- intersect(intedges,which.edge(tree,intaxa))
  #outtaxaedges <- intersect(intedges,which.edge(tree,outtaxa))
  #splitedge <- setdiff(intedges,union(intaxaedges,outtaxaedges))
  splitedge <- which.edge(tree,intaxa)
  
  if (length(intaxa)!=1){
    intaxaedges <- which.edge(tree,intaxa)
    outtaxaedges <- which.edge(tree,outtaxa)
    splitedge <- setdiff(seq(1,length(tree$edge[,1])),union(intaxaedges,outtaxaedges))
  }
  
  edgelabels(supportsites[i,2],splitedge,bg=mycolors[1],adj=c(.5,-1.1))
  edgelabels(supportsites[i,3],splitedge,bg=mycolors[2],adj=c(.5,.5))
  edgelabels(supportsites[i,4],splitedge,bg=mycolors[3],adj=c(.5,2.1))
  edgelabels(supportsites[i,5],splitedge,bg=mycolors[4],adj=c(.5,3.7))
}

#These labelling / legend things need adjustment sometimes
tiplabels()
#samplenumbers <- sort(sapply(tree$tip.label, FUN = function(x) strtoi(regmatches(x,regexpr('[0-9]+',x)))))
# samplenumbers <- sapply(samplenumbers, FUN = function(x) sub('M2c','M3a',x))
#tiplabels(sapply(tree$tip.label, FUN=function(x) substr(x,1,nchar(x)-1)),frame="none",adj=-.25)
#tiplabels(tree$tip.label,frame="none",adj=-.25)
#tiplabels(sapply(tree$tip.label, FUN=function(x) samplenumbers[x]),frame="none",adj=-.25)
legend('bottomleft',c('% of supporting sites','% of sites indicating a homozygous -> heterozygous \n mutation occured here','% of sites indicating a heterozygous -> homozygous \n mutation occured here','% of sites with a homozygote \n and a heterozygote on both sides'), fill = mycolors, cex = .9, y.intersp=2, bty='n')

dev.off()

cat('\ntotals\n======\n')
cat(paste0('support ',sum(supportsites[,2]),'\n'))
cat(paste0('heterize ',sum(supportsites[,3]),'\n'))
cat(paste0('homerize ',sum(supportsites[,4]),'\n'))
cat(paste0('incompatible ',sum(supportsites[,5]),'\n'))

q()