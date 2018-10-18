

# variant file for final set of variants
dng_variants = "/Users/roblanfear/Documents/github/somatic-variation/results/dng_variants.tsv"

# output file
outdir = "~/Documents/github/somatic-variation/variant_analyses/analysis3_physical_distribution/"
dir.create(outdir)


# function to return intverals between mutations
# mutations is a list of locations, where this:
# Chr1 (10bp): 3, 7
# Chr2 (10bp): 2, 6, 8
# is concatenated to this:
# c(3, 7, 12, 16, 18)
# and outputs this:
# 4, 4, 2
get_intervals = function(muts, chrlens){
    # muts: a list of mutation locations as above
    # chrlens: a list of chromosome lengths
    
    chrs = cumsum(chrlens)
    
    intervals = c()
    start = 0

    for(stop in chrs){
        print(paste(start, stop))
        m = muts[muts > start]
        m = m[m <= stop]
        m = sort(m)
        print(m)
        intervals = c(intervals, diff(m))
        print(intervals)
        start = stop
    }
    return(intervals)
}

simulate_mutations = function(chrlens, N=100, mindist = 500){
    
    #chrlens: a list of chromosome lengths
    #N: number of mutations to simulate
    #mindist: minimum distance one mutation has to be from another
    
    len = sum(chrlens)
    muts = c()
    while(length(muts)<N){
        
        new = sample(1:len, 1)
        if(length(muts)>0){ 
            if(min(abs(muts - new)) > mindist){
                muts = c(muts, new)
            }
        }else{
            muts = c(muts, new)
        }
    }
    
    print(muts)
    return(muts)
}


simulate_intervals = function(chrlens, N=100, mindist = 500){
    
    muts = simulate_mutations(chrlens, N = N, mindist = mindist)
    ints = get_intervals(muts, chrlens)

    return(ints)    
}


# chromosome lengths
chrlens = c(100, 100, 100, 100, 100, 100)
ints = simulate_intervals(chrlens, N = 100, mindist = 10)

# 100 replicates
ints = replicate(n = 100, simulate_intervals(chrlens, N = 100, mindist = 10))


# observed
d = read.delim(dng_variants)
d = d[,1:2]


# from the Egrandis reference
# cat Egra_chr1_chr11.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
#Chr01	44965197
#Chr02	59529170
#Chr03	83952244
#Chr04	41160059
#Chr05	76243064
#Chr06	57472304
#Chr07	54830020
#Chr08	72515979
#Chr09	39307835
#Chr10	37777128
#Chr11	44836791
chrlens = c(44965197, 59529170, 83952244, 41160059, 76243064, 57472304, 54830020, 72515979, 39307835, 37777128, 44836791)
chrsum = cumsum(chrlens)
d$mut = d$pos
chr = 2
for(l in chrsum){
    
    orig = which(d$chr == paste("scaffold_", chr, sep=""))
    
    print(l)
    
    d$mut[orig] = d$mut[orig] + l
    
    chr = chr+1
    
}

obs_ints = get_intervals(d$mut, chrlens)

rep_ints = ints = replicate(n = 100, simulate_intervals(chrlens, N = length(d$mut), mindist = 500))

mean_obs = mean(obs_ints)

repmean = mean(unlist(rep_ints))

mean_reps = as.data.frame(unlist(lapply(rep_ints, mean)))
names(mean_reps) = c('mean_interval', 'replicate')

p1 = ggplot(mean_reps, aes(x=mean_interval)) + geom_histogram(bins = 6) + 
    geom_vline(xintercept = mean_obs, colour = 'red', linetype = 'dashed', size = 2) + 
    geom_vline(xintercept = repmean, colour = 'red', size = 1)
    
pdf(file.path(outdir, "physical_distribution_figure.pdf"), width=7, height=5)
p1
dev.off()

rank = length(which(mean_reps$mean_interval < mean_obs))
p = (100-rank)/50 # two-tailed test

p

