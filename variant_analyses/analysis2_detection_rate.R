library(ggplot2)
library(stringr)

# input files
# must exist prior to running this script
callability = "/Users/roblanfear/Documents/github/somatic-variation/results/callability.txt"

outdir = "~/Documents/github/somatic-variation/variant_analyses/analysis2_detection_rate/"
dir.create(outdir)


####### Analysis 4: estimate false negative rate from in-silico mutations #########

d = read.delim(callability)

# sanity check - there should be 14 branches total
unique(d$samples_mutated)

# name the branches. We use Ancestor-Descendent notation
d$branch = NA

d$branch[which(d$samples_mutated == "M1")] = 'D->M1'
d$branch[which(d$samples_mutated == "M2")] = 'D->M2'
d$branch[which(d$samples_mutated == "M3")] = 'C->M3'
d$branch[which(d$samples_mutated == "M4")] = 'B->M4'
d$branch[which(d$samples_mutated == "M5")] = 'E->M5'
d$branch[which(d$samples_mutated == "M6")] = 'G->M6'
d$branch[which(d$samples_mutated == "M7")] = 'G->M7'
d$branch[which(d$samples_mutated == "M8")] = 'F->M8'

d$branch[which(d$samples_mutated == "M1,M2")] = 'C->D'
d$branch[which(d$samples_mutated == "M1,M2,M3")] = 'B->C'
d$branch[which(d$samples_mutated == "M1,M2,M3,M4")] = 'A->B'
d$branch[which(d$samples_mutated == "M6,M7")] = 'F->G'
d$branch[which(d$samples_mutated == "M6,M7,M8")] = 'E->F'
d$branch[which(d$samples_mutated == "M5,M6,M7,M8")] = 'A->E'

# sanity check
unique(d$branch)

pdf(file.path(outdir, "detection_rate_figure.pdf"))
ggplot(d, aes(x=depth)) + 
    geom_histogram(bins=100) + 
    facet_wrap(~mutation_recovered, ncol=1) +
    ggtitle("Detection rate", subtitle = "Coverage per site, with facets for whether or not an in-silico mutation was recovered")

ggplot(d, aes(x=depth)) + 
    geom_histogram(bins=100) + 
    facet_wrap(~mutation_recovered, ncol=1) + 
    xlim(0, 1000) +
    ggtitle("Detection rate, zoomed to x[0, 1000]", subtitle = "Coverage per site, with facets for whether or not an in-silico mutation was recovered")

ggplot(d, aes(x=depth)) + 
    geom_density(aes(colour = mutation_recovered)) + 
    facet_wrap(~branch) +
    xlim(0, 1000) +
    ggtitle("Detection rate per branch x[0, 1000]", subtitle = "Coverage per site, with facets for whether or not an in-silico mutation was recovered")

ggplot(d, aes(x=depth)) + 
    geom_density(aes(colour = mutation_recovered)) + 
    facet_wrap(~scaffold) +
    xlim(0, 1000) +
    ggtitle("Detection rate per scaffold x[0, 1000]", subtitle = "Coverage per site, with facets for whether or not an in-silico mutation was recovered")
    
dev.off()

# recovery rate per branch
t = table(d$mutation_recovered, d$branch)

t = data.frame(branch = colnames(t), recovered = t[2,], not.recovered = t[1,])
t$recovery.rate = t$recovered / (t$recovered + t$not.recovered)

t$branch.type = str_count(as.character(t$branch), "M")
t$branch.type[which(t$branch.type==1)] = "tip"
t$branch.type[which(t$branch.type==0)] = "internal"

# output results
write.table(t, file.path(outdir, "detection_rate_per_branch.tsv"), quote=FALSE, sep='\t', row.names = FALSE)

# save the overall recovery rate
overall.rec.rate = sum(t$recovered) / (sum(t$recovered) + sum(t$not.recovered))
fileConn<-file(file.path(outdir, "detection_rate_overall.txt"))
writeLines(c("The overall recovery rate from the in-silico mutations is:", as.character(overall.rec.rate)), fileConn)
close(fileConn)