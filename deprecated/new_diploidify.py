#!/usr/bin/env python3

"""
Take IUPAC alignment and remove sites that require 2 mutations (ie C/C -> T/T)
Also adds a variable site to reflect a mutation (ie R -> A/G and A -> A/A)
This is useful for using haploid tools with diploid data
Keep in mind we ignore sites with more than 2 alleles!
"""

import sys
import argparse
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment
import vcf

iupac = {
    'A' : 'AA',
    'C' : 'CC',
    'T' : 'TT',
    'G' : 'GG',
    'R' : 'AG',
    'Y' : 'CT',
    'S' : 'CG',
    'W' : 'AT',
    'K' : 'GT',
    'M' : 'AC',
    'N' : 'NN',
    '.' : '..',
    '-' : '--'
}

gaps = {
    'N' : 'NN',
    '.' : '..',
    '-' : '--'
}

def argparser():
    parser = argparse.ArgumentParser(description=
        """
        Create explicit alignment from a diploid sequence that utilizes IUPAC codes.
        Will remove sequences that are identical to previous sequences in the alignment.
        Will also remove any site that is not biallelic.
        Use "-t vcf" to perform these checks and remove nonvariable sites from a vcf.
        In VCF mode, this program can output a vcf or alignment.
        """)
    parser.add_argument("-i","--input", help="input file (default: stdin)")
    parser.add_argument("-t","--type", help="filetype (default: fasta). Any Biopython-compatible alignment format or \"vcf\".", default="fasta")
    parser.add_argument("-o","--output",help="output file (default: stdout)")
    parser.add_argument("-p","--outtype",help="output filetype (default: same as --type)")
    parser.add_argument("-v","--variable",help="only output variable sites",action="store_true")
    parser.add_argument("-s","--skip",help="keep single-letter notation but keep other checks",action="store_true")
    args = parser.parse_args()
    if not args.input: args.input = sys.stdin
    if not args.output: args.output = sys.stdout
    if not args.outtype: args.outtype = args.type
    return args

def remove_duplicate_seqs(alignment):
    """Remove duplicate sequences in the alignment; the replicates will all have identical sequences"""
    newseqs = list()
    returnseqs = list()
    for record in alignment:
        sequence = record.seq
        if sequence not in newseqs:
            newseqs.append(record.seq)
            # record.seq = Seq(str(record.seq).replace('N','-'), record.seq.alphabet)
            returnseqs.append(record)
    return MultipleSeqAlignment(returnseqs)

def is_biallelic(baselist):
    """Returns TRUE if list of bases is biallelic; false otherwise"""
    return unique_list_size(baselist) < 3

def filter_alignment(args):
    filein = args.input
    fileout = args.output
    filetype = args.type
    outtype = args.outtype
    variablesites = args.variable
    skip = args.skip

    seqs = AlignIO.read(filein,filetype)
    newalignment= [''] * len(seqs)
    for i in range(0,seqs.get_alignment_length()):
        baselist = list(seqs[:,i])
        if not is_biallelic(baselist):
            continue
        gts = [set(iupac[b.upper()]) for b in baselist] #each gt is a set, containing 1 or 2 members. this is a list of sets.
        c = gts[0].intersection(gts[1:]) #c is the common allele. this will be empty if there were 2 mutations (ie A/A -> T/T) so remove those samples
        if len(c) == 0:
            continue
        commonalleles = [c] * len(gts)
        variablealleles = [( (gt-c) if (gt-c) else gt ) for gt in gts] #if sample is homozygous for the common allele, its gt-c set will be empty, so its second allele is the same as its first

        if args.variable:
            commonalleles = [set()] * len(gts)
            if len(variablealleles[0].union(variablealleles[1:])) == 1:
                variablealleles = [set()] * len(gts)
            if commonalleles == variablealleles == [set()] * len(gts):
                continue

        combined = list()
        if skip:
            combined = baselist
        else:
            #convert to list of chars first
            commonalleles = [list(c)[0] if c else '' for c in commonalleles]
            variablealleles = [list(v)[0] if v else '' for v in variablealleles]
            combined = [c + v for c,v in zip(commonalleles, variablealleles)]

        newalignment = [newalignment[j] + combined[j] for j in range(0,len(combined))]

    newseqobjs = [SeqRecord(Seq(newalignment[l], IUPAC.unambiguous_dna), id=seqs[l].id, description='') for l in range(0,len(seqs))]
    newalnobj = MultipleSeqAlignment(newseqobjs)
    newalnobj = remove_duplicate_seqs(newalnobj)
    AlignIO.write(newalnobj,fileout,outtype)

def filter_vcf(args):
    fd = args.input
    of = args.output
    ot = args.outtype
    vcf_reader = (vcf.Reader(fsock=fd) if fd == sys.stdin else vcf.Reader(filename=fd))
    if ot == "vcf":
        vcf_writer = (vcf.Writer(of,vcf_reader) if of == sys.stdout else vcf.Writer(open(of,'w'),vcf_reader))
    else:
        newalignment = None
        samplenames = None

    for record in vcf_reader:
        gts = list()
        for sample in record.samples:
            sepchar = ('/' if not sample.phased else '|') #replacing this is not strictly necessary, but makes sense to do
            bases = sample.gt_bases
            if bases != None:
                gts.append(set(bases.replace(sepchar,''))) #sort to ensure A/T is not different from T/A
            else:
                gts.append(None)
                break
        if gts[-1] is None: continue
        c = gts[0].intersection(gts[1:]) #c is the common allele. this will be empty if there were 2 mutations (ie A/A -> T/T) so remove those samples
        if len(c) == 0:
            continue
        commonalleles = [c] * len(gts)
        variablealleles = [( (gt-c) if (gt-c) else gt ) for gt in gts] #if sample is homozygous for the common allele, its gt-c set will be empty, so its second allele is the same as its first

        if args.variable:
            commonalleles = [set()] * len(gts)
            if len(variablealleles[0].union(variablealleles[1:])) == 1:
                variablealleles = [set()] * len(gts)
            if commonalleles == variablealleles == [set()] * len(gts):
                continue

        combined = list()
        if skip:
            combined = baselist
        else:
            #convert to list of chars first
            commonalleles = [list(c)[0] if c else '' for c in commonalleles]
            variablealleles = [list(v)[0] if v else '' for v in variablealleles]
            combined = [c + v for c,v in zip(commonalleles, variablealleles)]

        if commonalleles == variablealleles == [set()] * len(gts):
                continue

        if ot == "vcf":
            vcf_writer.write_record(record)
        else:
            if not newalignment:
                newalignment = [''] * len(record.samples)
                samplenames = [x.sample for x in record.samples]
            newalignment = [newalignment[j] + combined[j] for j in range(0,len(combined))]

    if ot == "vcf":
        vcf_writer.close()
    else:
        newseqobjs = [SeqRecord(Seq(newalignment[l], IUPAC.unambiguous_dna), id=samplenames[l], description='') for l in range(0,len(newalignment))]
        newalnobj = MultipleSeqAlignment(newseqobjs)
        newalnobj = remove_duplicate_seqs(newalnobj)
        AlignIO.write(newalnobj,of,ot)

def main():
    args = argparser()
    if args.type != "vcf":
        filter_alignment(args)
    else:
        filter_vcf(args)

if __name__ == '__main__':
    main()



