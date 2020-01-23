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
import collections
import logging

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
    Given a vcf record and sample, return the sample genotype as a Counter.
    """
    return collections.Counter(record.genotype(sample).gt_alleles)

def test_vcf_record(record, replicates, n, m):
    """
    Given a vcf record and a dictionary of sample names -> list of
    replicate names, return a modified version of the record  if at
    least n replicates of every sample have matching genotypes. All
    samples will be modified to the most common genotype if it is >= n and
    >1. If n = 1, the sites will not be modified and the record will be
    returned. Returns None when the site fails the filter.
    """
    #if record.num_unknown > 0:
    #    logging.debug("Site with uncalled samples skipped")
    #    return None
    for s in replicates.keys():
        reps = replicates[s]
        gts = [tuple(get_sample_gt(record, r).items()) for r in reps]
        gtcounts = collections.Counter(gts)
        common_gt, common_num = gtcounts.most_common(1)[0]
        if m:
            possible_agreement = len(reps) - gtcounts[None]
            if possible_agreement > 2:
                n = possible_agreement/2
            else:
                n = possible_agreement
        if common_num < n:
            return None
        elif common_num == len(reps):
            continue
        elif n == 1 and not m:
            continue
        else:
            common_gt = list(collections.Counter(dict(common_gt)).elements())
            for r in reps:
                call = record.genotype(r)
                sepchar = call.gt_phase_char()
                call.gt_alleles = common_gt
                gtstr = sepchar.join(common_gt) if not None in common_gt else sepchar.join(['.','.'])
                call.data = call.data._replace(GT=gtstr)
    else:
        return record

def main():
    parser = argparse.ArgumentParser(
        description = "Remove sites that don't match the tree topology." )
    parser.add_argument("-v","--vcf", required = True, help = "Input vcf file")
    parser.add_argument("-p","--ped", required = True, help = "Input ped file")
    exclusive = parser.add_mutually_exclusive_group(required = True)
    exclusive.add_argument("-n", type = int, help = "Number of replicates required to match to keep site.")
    exclusive.add_argument("-m", action = 'store_true', help = "Majority Rule")
    args = parser.parse_args()

    samples = ['M' + str(s + 1) for s in range(8)]

    treestr = parse_pedfile(args.ped)
    replicates = parse_tree(treestr, samples)

    vcfr = vcf.Reader(filename = args.vcf) if args.vcf is not '-' else vcf.Reader(sys.stdin)
    vcfw = vcf.Writer(sys.stdout, vcfr)
    for record in vcfr:
        newrec = test_vcf_record(record, replicates, args.n, args.m)
        if newrec is not None:
            vcfw.write_record(record)
            #for BED format:
            #print(record.CHROM, record.start, record.end, sep = "\t")

if __name__ == "__main__":
    main()
