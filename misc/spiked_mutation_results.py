#!/usr/bin/env python2
import pysam
import vcf

def load_muts():
    sample_files = ["mut_files/mut_M" + str(k) + ".txt" for k in range(1,9)] #get every file
    muts = []
    for fn in sample_files:
        with open(fn) as fh:
            d = [[l.split()[0], int(l.split()[1])] for l in fh] #a list of every [scaffold,site] in the sample
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
        if col.reference_pos != position - 1:
            continue # -1 because pysam are 0-based coordinates but samtools is 1-based
        else:
            return col.nsegments #just first position in pileup
    else:
        raise ValueError("This pileup is invalid!" + str(chr) + ':' + str(position))

"""
scaffold, site, original_genotype, mutated_genotype, depth, branch_mutated, samples_mutated, mutation_recovered
"""
def generate_table_line(line, muts, vcf, sam):
    l = line.rstrip().split()
    loc = [l[0], int(l[1])]
    gt = l[2:4]
    togt = set(l[4:6])
    depth = position_depth(sam, l[0], l[1])
    mutated_samples = ["M" + str(i) for i in range(1,9) if loc in muts[i-1]] #sample names are 1-based
    if mutated_samples == []:
        return None

    recovered = False
    if loc[0] in vcf:
        chrd = vcf[loc[0]]
        if loc[1] in chrd:
            record = chrd[loc[1]]
            for s in mutated_samples:
                vcfgt = record.genotype(s + "a").gt_bases # add a because the samples are M1a, M2a, etc.
                if vcfgt == None:
                    break
                vcfg = set(vcfgt.split("/"))
                if vcfg != togt:
                    break
            else: #if we make it through the loop, we recovered the mutation
                recovered = True

    return loc[0], str(loc[1]), ''.join(gt), ''.join(togt), str(depth), ','.join(mutated_samples), str(recovered)

def main():
    muts = load_muts() #open all mut files and load into list
    vcf = load_vcf()
    sam = load_sam()

    mutfile = "mut_files/mutfile.txt" #iterate through full list and generate table
    print "\t".join(["scaffold","site","original_genotype","mutated_genotype","depth","samples_mutated", "mutation_recovered"])
    with open(mutfile) as fh:
        for l in fh:
            outline = generate_table_line(l, muts, vcf, sam)
            if outline != None:
                print "\t".join(outline)

if __name__ == '__main__':
    main()



