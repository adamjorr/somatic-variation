#!/usr/bin/env python3
#ensure the "mut_files" directory exists in the working directory.
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

# def load_bed(bedfile):
#     bed = BedTool(bedfile)
#     bedranges = dict() #dictionary with bedranges[chr] = [range(x,y), range(z,a), ...]
#     for interval in bed:
#         l = bedranges.setdefault(interval.chrom, []) #list
#         bedranges[interval.chrom] = l + [range(interval.start+1, interval.stop+1)] #switch BED to 1-based
#     return bedranges

def load_repeats(repeatfile):
    with open(repeatfile) as fh:
        repeats = [[l.split()[0], int(l.split()[1])] for l in fh]
    return repeats

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
def generate_table_line(line, muts, vcf, sam, repeats, dng = False, norepfilter = False):
    l = line.rstrip().split()
    loc = [l[0], int(l[1])]
    gt = l[2:4]
    togt = set(l[4:6])
    depth = position_depth(sam, l[0], l[1])
    #in dng, there is SM/M1, SM/M2, etc. Otherwise, the samples are M1a, M2a, M3a ... etc
    if norepfilter:
        mutated_samples = ["M" + str(i) + j for i in range(1,9) for j in ["a","b","c"] if loc in muts[i - 1]]
    elif dng:
        mutated_samples = ["SM/M" + str(i) for i in range(1,9) if loc in muts[i-1]] #sample names are 1-based
    else:
        mutated_samples = ["M" + str(i) + "a" for i in range(1,9) if loc in muts[i-1]]
    if mutated_samples == []:
        return None

    inrepeat = (loc in repeats)

    recovered = False
    nunknown = 0
    if loc[0] in vcf:
        chrd = vcf[loc[0]]
        if loc[1] in chrd:
            record = chrd[loc[1]]
            nunknown = record.num_unknown
            # gts = [record.genotype(s).gt_bases for s in mutated_samples]
            # if not None in gts:
            #     gts = [set(g.split("/")) for g in gts]
            if norepfilter:
                #no rep filter
                vcfgts = [record.genotype(s).gt_bases for s in mutated_samples]
                vcfgts = list(zip(*[iter(vcfgts)]*3)) #list of tuples of len 3
                for s in vcfgts:
                    if all([x == None for x in s]):
                        break
                    if all([set(x.split("/")) != togt for x in s]):
                        break
                else: #made it thru the loop, count it as recovered
                    recovered = True
            else:
                for s in mutated_samples:
                    vcfgt = record.genotype(s).gt_bases
                    if vcfgt == None:
                        break
                    vcfg = set(vcfgt.split("/"))
                    if vcfg != togt:
                        break
                else: #if we make it through the loop, we recovered the mutation
                    recovered = True

    return loc[0], str(loc[1]), ''.join(gt), ''.join(togt), str(depth), ','.join(mutated_samples), str(recovered), str(inrepeat), str(nunknown)

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
    # parser.add_argument("-b","--bedfile", help = "BED file containing repeat regions")
    parser.add_argument("-r","--repeatfile", help = "File containing mutated locations in repeat regions.")
    parser.add_argument("--dng", action = 'store_true', help = "Expect dng-style VCF, where the sample is coded as SM/MX rather than MXa")
    parser.add_argument("--norepfilter", action = 'store_true', help = "The data is not replicate filtered. If any replicate is correct, the mutation will count as detected.")

    args = parser.parse_args()
    return args

def main():
    args = argparser()
    vcffile = args.vcffile
    samfile = args.samfile
    bedfile = args.repeatfile
    muts = load_muts() #open all mut files and load into list
    vcf = load_vcf(vcffile) if vcffile is not None else load_vcf()
    sam = load_sam(samfile) if samfile is not None else load_sam()
    repeats = load_repeats(bedfile)

    mutfile = "mut_files/mutfile.txt" #iterate through full list and generate table
    print("\t".join(["scaffold","site","original_genotype","mutated_genotype","depth","samples_mutated", "mutation_recovered", "in_repeat_region", "num_unknown_gts"]))
    with open(mutfile) as fh:
        for l in fh:
            outline = generate_table_line(l, muts, vcf, sam, repeats, args.dng, args.norepfilter)
            if outline != None:
                print("\t".join(outline))

if __name__ == '__main__':
    main()



