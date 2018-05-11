library(phangorn)
library(ggtree)
library(ggplot2)

# input files
# must exist prior to running this script
alignment = "~/Documents/github/somatic-variation/variant_analyses/analysis0_alignments/aln_0gaps.fa"
true.tree = "~/Documents/github/somatic-variation/variant_analyses/tree_em1_physical_structure.phy"

# make output directory
outdir = "~/Documents/github/somatic-variation/variant_analyses/analysis4_mutation_spectrum/"
dir.create(outdir)


################ Analysis 4 #####################
#       calculate the mutation spectrum         #

# load the alignment
EM1 = read.phyDat(alignment, format="fasta", type="DNA")

# load the physical tree, and calculate MP branch lengths
tree = read.tree(true.tree)
tree = unroot(tree)
tree = acctran(tree, EM1)

# sanity check: make sure this is the expected tree and labelling scheme
plot(tree)
nodelabels(text = tree$node.label)
ggtree(tree) + geom_label(aes(label=label)) #prettier!


# Now we identify all of the mutations identified by the MP analysis
# and build those into a detailed dataframe so we can look at all mutations
anc.p = ancestral.pars(tree, EM1, "ACCTRAN")
anc.p = phangorn:::ancestral2phyDat(anc.p)
charMatrix = as.character(anc.p)
dim(charMatrix)
row.names(charMatrix)

res = NULL
for(i in 1:nrow(tree$edge)){
  e1 <- tree$edge[i, 1]
  e2 <- tree$edge[i, 2]
  ind <- which(charMatrix[e1,] != charMatrix[e2,])
  if(length(ind)>0)res = rbind(res, cbind(e1,e2,charMatrix[e1,ind], charMatrix[e2,ind], ind))
}

res = as.data.frame(res)

names(res) = c("anc.node", "des.node", "anc.base", "des.base", "aln.pos")

# now we make a column of mutation names
# accounting for complementarity in the standard way
res$mut = paste(res[,3], res[,4])
res$mutrev = res$mut

res$mutrev[which(res$mut=='a g')] = 'at->gc'
res$mutrev[which(res$mut=='t c')] = 'at->gc'

res$mutrev[which(res$mut=='g a')] = 'gc->at'
res$mutrev[which(res$mut=='c t')] = 'gc->at'

res$mutrev[which(res$mut=='t g')] = 'ta->gc'
res$mutrev[which(res$mut=='a c')] = 'ta->gc'

res$mutrev[which(res$mut=='a t')] = 'at->ta'
res$mutrev[which(res$mut=='t a')] = 'at->ta'

res$mutrev[which(res$mut=='g t')] = 'gc->ta'
res$mutrev[which(res$mut=='c a')] = 'gc->ta'

res$mutrev[which(res$mut=='g c')] = 'gc->cg'
res$mutrev[which(res$mut=='c g')] = 'gc->cg'

# order the factor levels for the plot
res$mutrev = factor(res$mutrev, levels = c("at->gc", "gc->at", "ta->gc", "at->ta", "gc->ta", "gc->cg"))


# now we do some long-winded renaming of things to 
# annotate exaclty which mutations happened on which 
# branches in the tree
# the renaming scheme can be seen by examining these two trees
# tree1
plot(tree, show.tip.label = F)
nodelabels()
tiplabels()

# tree2
plot(tree, show.tip.label = T)
nodelabels(text = tree$node.label)

# set up dummy variables
res$anc.node = as.integer(as.character(res$anc.node))
res$des.node = as.integer(as.character(res$des.node))
res$anc.nodelab = res$anc.node
res$des.nodelab = res$des.node

# label the branches properly
res$anc.nodelab[which(res$anc.node==1)] = 'M1'
res$des.nodelab[which(res$des.node==1)] = 'M1'
res$anc.nodelab[which(res$anc.node==2)] = 'M2'
res$des.nodelab[which(res$des.node==2)] = 'M2'
res$anc.nodelab[which(res$anc.node==3)] = 'M3'
res$des.nodelab[which(res$des.node==3)] = 'M3'
res$anc.nodelab[which(res$anc.node==4)] = 'M4'
res$des.nodelab[which(res$des.node==4)] = 'M4'
res$anc.nodelab[which(res$anc.node==8)] = 'M5'
res$des.nodelab[which(res$des.node==8)] = 'M5'
res$anc.nodelab[which(res$anc.node==6)] = 'M6'
res$des.nodelab[which(res$des.node==6)] = 'M6'
res$anc.nodelab[which(res$anc.node==5)] = 'M7'
res$des.nodelab[which(res$des.node==5)] = 'M7'
res$anc.nodelab[which(res$anc.node==7)] = 'M8'
res$des.nodelab[which(res$des.node==7)] = 'M8'


res$anc.nodelab[which(res$anc.node==12)] = 'E'
res$des.nodelab[which(res$des.node==12)] = 'E'
res$anc.nodelab[which(res$anc.node==11)] = 'D'
res$des.nodelab[which(res$des.node==11)] = 'D'
res$anc.nodelab[which(res$anc.node==10)] = 'C'
res$des.nodelab[which(res$des.node==10)] = 'C'
res$anc.nodelab[which(res$anc.node==9)] = 'B'
res$des.nodelab[which(res$des.node==9)] = 'B'
res$anc.nodelab[which(res$anc.node==13)] = 'F'
res$des.nodelab[which(res$des.node==13)] = 'F'
res$anc.nodelab[which(res$anc.node==14)] = 'G'
res$des.nodelab[which(res$des.node==14)] = 'G'

res$branch = paste(res$anc.nodelab, res$des.nodelab, sep="->")

# mutation spectrum for the whole dataset
ggplot(res, aes(x=mutrev)) + geom_histogram(stat="count") +
    xlab("mutation type")

# transition/transversion ratio
t = table(res$mutrev)
ts = t["at->gc"] + t["gc->at"]
ts.tv =  ts / ( sum(t) - ts ) 
ts.tv
