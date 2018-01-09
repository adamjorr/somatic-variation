#!/usr/bin/env python2
from pybedtools import BedTool
import pybedtools
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio import SeqIO
import argparse
import sys
import string

def alt_AAs(site, codon, tab):
    """
    site = one of (0,1,2)
    """
    table = unambiguous_dna_by_id[tab]
    altcodons = [ codon[:site] + k + codon[site+1:] for k in ['A','T','C','G'] ]
    return [table.forward_table[c] if c in table.forward_table else '*' for c in altcodons]

def is_degenerate(site, codon, table):
    return (len(set(alt_AAs(site,codon,table))) == 1)

def degenerate_pos(seq, table = 1):
    """Get 0-based positions consisting of the x-fold degenerated codons only."""
    data = []
    for i in range(0, len(seq), 3):
        codon = str(seq[i:i + 3])
        if (len(codon) == 3):
            for j in range(3):
                if is_degenerate(j,codon,table):
                    data.append(i + j)
    return data

def argparser():
    parser = argparse.ArgumentParser(description=
        """
        Find 4-fold degenerate sites given a reference and BED file.
        """)
    parser.add_argument("-i","--input", help="input BED file (default: stdin)")
    parser.add_argument("-o","--output",help="output file (default: stdout)")
    parser.add_argument("-r","--reference",help="reference file", required = True)
    args = parser.parse_args()
    if not args.input: args.input = sys.stdin
    if not args.output: args.output = sys.stdout
    return args

# def getref(args):
#     """
#     return dictionary of {chr : string}
#     """
#     return SeqIO.to_dict(SeqIO.parse(args.reference,"fasta"))

def getCDSs(bedfilename, reffilename, strand):
    """
    return iterator of coding sequences
    """
    bed = BedTool(bedfilename)
    bed = bed.filter(lambda x: x.strand == strand)
    fasta = reffilename
    bed = bed.sequence(fi=fasta, s=True)
    return SeqIO.parse(bed.seqfn, "fasta")

def findSites(bedfilename, reffilename, strand):
    """
    if reverse strand, absolute position is len(seq) - pos
    """
    sites = []
    for record in getCDSs(bedfilename, reffilename, strand):
        pos = degenerate_pos(record.seq)
        refchr, refpos = record.id.split(':', 1)
        refpos, _ = refpos.split('-', 1)
        refpos = int(refpos)
        for k in pos:
            p = (refpos + k if strand == "+" else refpos + len(record.seq) - k - 1)
            sites.append('\t'.join(refchr,str(p),str(p+1)))
    return sites

def main():
    args = argparser()
    forwardsites = findSites(args.input, args.reference, "+")
    revsites = findSites(args.input, args.reference, "-")
    try:
        with open(args.output,'wb') as of:
            for l in forwardsites + revsites:
                of.write(l + "\n")
    except TypeError:
        for l in forwardsites + revsites:
            args.output.write(l + "\n")
    pybedtools.cleanup(remove_all=True)


if __name__ == '__main__':
    main()
