#!/usr/bin/env python2
import pysam
import vcf
import sys
import argparse
import pybedtools
from pybedtools import BedTool

def load_muts():
    sample_files = ["mut_files/mut_M" + str(k) + ".txt" for k in range(1,9)] #get every file
    muts = []
    for fn in sample_files:
        with open(fn) as fh:
            d = [[l.split()[0], int(l.split()[1])] for l in fh] #a list of every [scaffold,site] in the sample
            muts.append(d)
    return muts #a list of lists of lists

def load_vcf(vcffile = "repfiltered_only_mutated.vcf.gz"):
    vcfreader = vcf.Reader(filename = vcffile)
    vcflocs = dict()
    for record in vcfreader:
        d = vcflocs.setdefault(record.CHROM, dict()) #dictionary with vcflocs{chr} = {pos : record}
        d[record.POS] = record
        vcflocs[record.CHROM] = d
    return vcflocs

def load_sam(samfile = "fixed_alignment.bam"):
    return pysam.AlignmentFile(samfile, "rb")

def load_bed(bedfile):
    bed = BedTool(bedfile)
    bedranges = dict() #dictionary with bedranges[chr] = [range(x,y), range(z,a), ...]
    for interval in bed:
        l = bedranges.setdefault(interval.chrom, []) #list
        bedranges[interval.chrom] = l + [range(interval.start+1, interval.stop+1)] #switch BED to 1-based
    return bedranges

def position_depth(samfile, chr, position):
    regionstr = str(chr) + ':' + str(position) + '-' + str(position)
    cols = [c for c in samfile.pileup(region = regionstr, truncate = True)]
    if len(cols) == 0:
        return 0
    elif len(cols) > 1:
        raise ValueError("There were too many columns returned for this position:" + regionstr)
    elif cols[0].reference_pos != int(position) - 1:
        raise ValueError("This pileup is for the wrong position:" + regionstr)
    else:
        return cols[0].nsegments

"""
scaffold, site, original_genotype, mutated_genotype, depth, branch_mutated, samples_mutated, mutation_recovered
"""
def generate_table_line(line, muts, vcf, sam, bed):
    l = line.rstrip().split()
    loc = [l[0], int(l[1])]
    gt = l[2:4]
    togt = set(l[4:6])
    depth = position_depth(sam, l[0], l[1])
    mutated_samples = ["M" + str(i) for i in range(1,9) if loc in muts[i-1]] #sample names are 1-based
    if mutated_samples == []:
        return None

    inrepeat = False
    if loc[0] in bed:
        inrepeat = any([loc[1] in r for r in bed[loc[0]]])

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

    return loc[0], str(loc[1]), ''.join(gt), ''.join(togt), str(depth), ','.join(mutated_samples), str(recovered), str(inrepeat)

def argparser():
    parser = argparse.ArgumentParser(description =
        """
        Generate a table to estimate callability using a mutation table, VCF and SAM file.
        The mutation tables must be in a subdirectory called mut_files.
        The output table is printed to STDOUT.
        """)
    #parser.add_argument("-m","--mutfile")
    parser.add_argument("-v","--vcffile", help = "vcf file")
    parser.add_argument("-s","--samfile", help = "sam/bam file")
    parser.add_argument("-b","--bedfile", help = "BED file containing repeat regions", )

    args = parser.parse_args()
    return args

def main():
    args = argparser()
    vcffile = args.vcffile
    samfile = args.samfile
    bedfile = args.bedfile
    muts = load_muts() #open all mut files and load into list
    vcf = load_vcf(vcffile) if vcffile is not None else load_vcf()
    sam = load_sam(samfile) if samfile is not None else load_sam()
    bed = load_bed(bedfile)

    mutfile = "mut_files/mutfile.txt" #iterate through full list and generate table
    print "\t".join(["scaffold","site","original_genotype","mutated_genotype","depth","samples_mutated", "mutation_recovered"])
    with open(mutfile) as fh:
        for l in fh:
            outline = generate_table_line(l, muts, vcf, sam)
            if outline != None:
                print "\t".join(outline)

if __name__ == '__main__':
    main()



