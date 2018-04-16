#!/usr/bin/env python3

import vcf
import pysam
import primer3

def get_ref_dict(reffilename):
    reffile = pysam.FastaFile(reffilename)
    return {r : reffile.fetch(region=r) for r in reffile.references}

def make_primers_in_vcf(vcffile, referencedict):
    vcfreader = vcf.Reader(filename = vcffile)
    for record in vcfreader:
        print(make_primer_for_record(record, referencedict))
        quit()
    return vcflocs

def make_primer_for_record(record, referencedict):
    chrom = record.CHROM
    pos = record.POS
    seqid = str(chrom) + str(pos)
    seq_args = {
        'SEQUENCE_ID' : seqid,
        'SEQUENCE_TEMPLATE' : referencedict[chrom],
        'SEQUENCE_TARGET' : [pos, 1]
    }
    # global_args = {
    #     'PRIMER_TASK' : 'pick_discriminative_primers',
    #     'PRIMER_FIRST_BASE_INDEX' : 1,
    #     'PRIMER_PICK_ANYWAY' : 1,
    #     'PRIMER_PRODUCT_SIZE_RANGE' : [[36,300]]
    # }
    global_args = {
        'PRIMER_TASK' : 'pick_sequencing_primers',
        'PRIMER_FIRST_BASE_INDEX' : 1,
        'PRIMER_PICK_ANYWAY' : 1
    }

    # libdict = {k:v for k,v in referencedict.items() if k is not chrom}
    return primer3.bindings.designPrimers(seq_args = seq_args, global_args = global_args)

def main():
    refdict = get_ref_dict('/storage/2017-03-15-eucalyptus/e_mel_3/e_mel_3.fa')
    make_primers_in_vcf('/storage/2017-03-15-eucalyptus/reeds_calls/filtered_calls/adams_calls.vcf.gz',refdict)

if __name__ == '__main__':
    main()
