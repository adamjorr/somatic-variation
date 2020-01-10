#!/usr/bin/env python3
"""
This script reads a PED file and a VCF and extracts the sites that have
the same genotype in N replicates as specified in the PED file.

The output is a VCF.
"""

import vcf
import Bio
import Bio.Phylo
import argparse
import sys

def parse_pedfile(filename):
    """
    The pedfile starts with a header line. The tree is in the 5th column.

    Opens filename and returns the tree as a string.
    """
    with open(filename) as fh:
        fh.readline() #skip header
        return fh.readline().split()[-1]

def parse_tree(treestr, samples):
    """
    Parse a newick tree. Our trees all have labels and branch lengths.

    Samples is the list of sample names to look for.

    Returns a dictionary of sample names -> list of replicate names.
    """
    replicates = {}
    tree = next(Bio.Phylo.NewickIO.Parser.from_string(treestr).parse())
    for c in tree.find_clades():
        if c.name and c.name in samples:
            replicates[c.name] = [d.name for d in c.clades]
    return replicates

def get_sample_gt(record, sample):
    """
    Given a vcf record and sample, return the sample genotype as a frozenset.
    """
    return frozenset(record.genotype(sample).gt_alleles)

def test_vcf_record(record, replicates, n):
    """
    Given a vcf record and a dictionary of sample names -> list of
    replicate names, return True if at least n replicates of every sample
    have matching genotypes. Otherwise return False.
    """
    if record.num_unknown > 0:
        return False
    for s in replicates.keys():
        reps = replicates[s]
        gts = [get_sample_gt(record, r) for r in reps]
        gtcounts = collections.Counter(gts)
        if gtcounts.most_common(1)[0][1] < n:
            return False
    else:
        return True

def main():
    parser = argparse.ArgumentParser(
        description = "Remove sites that don't match the tree topology." )
    parser.add_argument("-v","--vcf", required = True, help = "Input vcf file")
    parser.add_argument("-p","--ped", required = True, help = "Input ped file")
    parser.add_argument("-n", required = True, type = int, help = "Number of replicates required to match to keep site.")
    args = parser.parse_args()

    samples = ['M' + str(s + 1) for s in range(8)]
    replicates = parse_tree(parse_pedfile(args.ped), samples)

    vcfr = vcf.Reader(filename = args.vcf)
    vcfw = vcf.Writer(sys.stdout, vcfr)
    for record in vcfr:
        if test_vcf_record(record, replicates, args.n):
            vcfw.write_record(record)
            #for BED format:
            #print(record.CHROM, record.start, record.end, sep = "\t")

if __name__ == "__main__":
    main()
