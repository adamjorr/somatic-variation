import pybedtools
from pybedtools import BedTool
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio import SeqIO
import argparse

def alt_AAs(site, codon, tab):
    """
    site = one of (0,1,2)
    """
    table = unambiguous_dna_by_id[tab]
    allcodons = [codon[site]=k for k in ['a','t','c','g']]
    return [table[c] for c in altcodons]

def is_degenerate(site, codon, table):
    return (len(unique(alt_AAs(site,codon,table))) == 1)

def degenerate_pos(seq, table = "Standard"):
    """Get 1-based positions consisting of the x-fold degenerated codons only."""
    sites = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3].tostring()
        for j in range(3):
            if is_degenerate(j,codon,table):
                data.append(i + j + 1)
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
    bed = pybedtools.BedTool(bedfilename)
    bed = bed.filter(lambda x: x.strand == strand)
    fasta = reffilename
    bed = bed.sequence(fi=fasta, s=True)
    return SeqIO.parse(bed.seqfn, "fasta")

def findSites(bedfilename, reffilename, strand):
    """
    if reverse strand, absolute position is len(seq) + 1 - pos
    """
    sites = []
    for record in getCDSs(bedfilename, reffilename, strand):
        pos = degenerate_pos(record.seq)
        refchr, refpos = string.split(record.id, ':', 1)
        if strand == "+":
            sites.extend([refchr + ':' + (int(refpos)+k).tostring() for k in pos])
        else:
            sites.extend([refchr + ':' + (int(refpos)+(len(record.seq)+1 - k)) for k in pos])
    return sites

def main():
    args = argparser()
    forwardsites = findSites(args.input, args.reference, "+")
    revsites = findSites(args.input, args.reference, "-")
    with open(args.output,'wb') as of:
        for l in forwardsites + revsites:
            of.write(l + "\n")


if __name__ == '__main__':
    main()
