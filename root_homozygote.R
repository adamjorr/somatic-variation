#!/usr/bin/Rscript
#Rscript root_homozygote.R homozygosity_table.txt support_table.txt tree.nwk

library(ape)
library(phytools)
mycolors <- c('lightgreen','yellow','purple','red')
pdf('tree.pdf',8,8)

#args <- c('/home/minion/documents/eucalyptus/stats_on_tree/zygosity.txt','/home/minion/documents/eucalyptus/stats_on_tree/supp_table.txt','/home/minion/documents/eucalyptus/stats_on_tree/RAxML_bipartitions.diplified')
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
#tiplabels(sapply(tree$tip.label, FUN=function(x) substr(x,1,nchar(x)-1)),frame="none",adj=-.25)
tiplabels(tree$tip.label,frame="none",adj=-.25)
legend('topright',c('% of supporting sites','% of sites indicating a homozygous -> heterozygous \n mutation occured here','% of sites indicating a heterozygous -> homozygous \n mutation occured here','% of sites with a homozygote \n and a heterozygote on both sides'), fill = mycolors, cex = .9, y.intersp=2, bty='n')

dev.off()

cat('\ntotals\n======\n')
cat(paste0('support ',sum(supportsites[,2]),'\n'))
cat(paste0('heterize ',sum(supportsites[,3]),'\n'))
cat(paste0('homerize ',sum(supportsites[,4]),'\n'))
cat(paste0('incompatible ',sum(supportsites[,5]),'\n'))




q()
#This labels tree edges with the difference in number of homozygous sites.
for (i in seq(1,length(othersamples)-1)){
  intnode1 = tree$edge[which(tree$edge[,2]==largenodenum)]
  intnode2 = tree$edge[which(tree$edge[,1]==intnode1 & tree$edge[,2] != largenodenum),2]
  othernodenum = tree$edge[which(tree$edge[,1]==intnode2 & tree$edge[,2] <= length(tree$tip.label)),2]
  edgeidx = which(tree$edge[,1] == intnode1 & tree$edge[,2] == intnode2)
   
  smallname = tree$tip.label[largenodenum]
  if (length(othernodenum) == 2){
    othername1 = tree$tip.label[othernodenum[1]]
    othername2 = tree$tip.label[othernodenum[2]]
    otherval1 = counts[which(counts$V1 == othername1),][[2]]
    otherval2 = counts[which(counts$V1 == othername2),][[2]]
    if(otherval1 < otherval2){
      othernodenum = othernodenum[1]
    }
    else{
      othernodenum = othernodenum[2]
    }
  }
  othername = tree$tip.label[othernodenum]
  smallval = counts[which(counts$V1 == smallname),][[2]]
  otherval = counts[which(counts$V1 == othername),][[2]]
  dval = otherval - smallval
  edgelabels(dval,edgeidx)
  largenodenum = othernodenum
  
}
