#!/usr/bin/env python2
import pysam
import vcf

def load_muts():
    sample_files = ["mut_files/mut_M" + str(k) + ".txt" for k in range(1,9)] #get every file
    muts = []
    for fn in sample_files:
        with open(fn) as fh:
            d = [l.split()[0:2] for l in fh] #a list of every [scaffold,site] in the sample
            muts.append(d)
    return muts #a list of lists of lists

def load_vcf():
    vcffile = "repfiltered_only_mutated.vcf.gz"
    vcfreader = vcf.Reader(filename = vcffile)
    vcflocs = dict()
    for record in vcfreader:
        d = vcflocs.setdefault(record.CHROM, dict()) #dictionary with vcflocs{chr} = {pos : record}
        d[record.POS] = record
        vcflocs[record.CHROM] = d
    return vcflocs

def load_sam():
    samfile = "fixed_alignment.bam"
    return pysam.AlignmentFile(samfile, "rb")

def position_depth(samfile, chr, position):
    regionstr = str(chr) + ':' + str(position)
    for col in samfile.pileup(region = regionstr):
        return col.n #just first position in pileup

"""
scaffold, site, original_genotype, mutated_genotype, depth, branch_mutated, samples_mutated, mutation_recovered
"""
def generate_table_line(line, muts, vcf, sam):
    l = line.rstrip().split()
    loc = l[0:2]
    gt = l[2:4]
    togt = set(l[4:6])
    depth = position_depth(sam, l[0], l[1])
    mutated_samples = ["M" + str(i) for i in range(1,9) if loc in muts[i-1]] #sample names are 1-based
    
    recovered = False
    if loc[0] in vcf:
        chrd = vcf[loc[0]]
        try:
            if loc[1] in chrd:
                record = chrd[loc[1]]
                vcfgt = [set(record.genotype(s + "a").gt_bases.split("/")) for s in mutated_samples] # add a because the samples are M1a, M2a, etc.
                eq = [g == togt for g in vcfgt] #are all the genotypes equal to the togt?
                if False not in eq: #if so, we recovered the mutation
                    recovered = True
        except IndexError:
            print line;
            print l;
            print loc;
            raise

    return loc[0], loc[1], ''.join(gt), ''.join(togt), depth, ','.join(mutated_samples), recovered

def main():
    muts = load_muts() #open all mut files and load into list
    vcf = load_vcf()
    sam = load_sam()

    mutfile = "mut_files/mutfile.txt" #iterate through full list and generate table
    with open(mutfile) as fh:
        for l in fh:
            outline = generate_table_line(l, muts, vcf, sam)

if __name__ == '__main__':
    main()



