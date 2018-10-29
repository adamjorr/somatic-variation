#!/usr/bin/env python3

import vcf
import pysam
import primer3
import argparse
import sys

def get_ref_dict(reffilename):
    reffile = pysam.FastaFile(reffilename)
    return {r : reffile.fetch(region=r) for r in reffile.references}

def make_primers_in_vcf(vcfreader, referencedict, exclusions):
    print('\t'.join(['chr', 'pos', 'left', 'right', 'leftpos', 'rightpos']))
    for record in vcfreader:
        res = make_primer_for_record(record, referencedict, exclusions)
        res = short_result(res)
        if not res['left'] or not res['right']:
            print("Could not create a primer pair for coordinate", res['chr'], res['pos'], file = sys.stderr)
            continue
        line = array_result(res)
        print('\t'.join(line))

def make_primer_for_record(record, referencedict, exclusions):
    chrom = record.CHROM
    pos = record.POS
    seqid = str(chrom) + ":" + str(pos)

    seq_args = {
        'SEQUENCE_ID' : seqid,
        'SEQUENCE_TEMPLATE' : referencedict[chrom],
        'SEQUENCE_TARGET' : [pos, 1]
    }

    if exclusions is not None:
        ex = exclusions[chrom]
        seq_args['SEQUENCE_EXCLUDED_REGION'] = [e for e in ex if e is not [pos,1]]

    global_args = {
        'PRIMER_TASK' : 'pick_sequencing_primers',
        'PRIMER_FIRST_BASE_INDEX' : 1
    }

    # global_args = {
    #     'PRIMER_TASK' : 'pick_discriminative_primers',
    #     'PRIMER_FIRST_BASE_INDEX' : 1,
    #     'PRIMER_PICK_ANYWAY' : 1,
    #     'PRIMER_PRODUCT_SIZE_RANGE' : [[36,300]]
    # }
    # libdict = {k:v for k,v in referencedict.items() if k is not chrom}
    result = primer3.bindings.designPrimers(seq_args = seq_args, global_args = global_args)
    result['chr'] = chrom
    result['pos'] = pos
    return result

def short_result(p3dict):
    return {
        'chr' : p3dict.get('chr'),
        'pos' : p3dict.get('pos'),
        'left' : p3dict.get('PRIMER_LEFT_0_SEQUENCE',''),
        'right' : p3dict.get('PRIMER_RIGHT_0_SEQUENCE',''),
        'leftpos' : p3dict.get('PRIMER_LEFT_0',[''])[0],
        'rightpos' : p3dict.get('PRIMER_RIGHT_0',[''])[0]
        }

def array_result(res):
    return [str(res[k]) for k in ['chr', 'pos', 'left', 'right', 'leftpos', 'rightpos']]

def parse_args():
    parser = argparse.ArgumentParser(description = 
        """Create primers using primer3 from a given reference
        to sequence sites in a given vcf.
        """
        )
    parser.add_argument("-v", "--vcf", help = "VCF containing sites to amplify", default = sys.stdin)
    parser.add_argument("-r", "--reference", help = "Reference to create primers from")
    parser.add_argument("-e", "--exclude", help = """A vcf containing potentially variable sites to
                                                     exclude from overlapping the primers.""")
    args = parser.parse_args()
    return args

def load_files(vcffile, reffile):
    reader = vcf.Reader(filename = vcffile)
    if reffile is None:
        print("Reference not given. Attempting to find one in the VCF header.", file = sys.stderr)
        reffile = reader.metadata['reference']
        if reffile.startswith('file://'):
            reffile = reffile.replace('file://','',1)
        print("Using", reffile, "as reference.", file = sys.stderr )
    refdict = get_ref_dict(reffile)
    return reader, refdict

def load_exclusions(exclusionfile):
    reader = vcf.Reader(filename = exclusionfile)
    d = dict()
    for record in reader:
        d.setdefault(record.CHROM, list()).append([record.POS,1])
    return d

def main():
    args = parse_args()
    vcfreader, refdict = load_files(args.vcf, args.reference)
    if args.exclude is not None:
        exclusions = load_exclusions(args.exclude)
    else:
        exclusions = None
    make_primers_in_vcf(vcfreader,refdict,exclusions)

if __name__ == '__main__':
    main()
