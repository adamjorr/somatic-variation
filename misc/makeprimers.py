#!/usr/bin/env python3

import vcf
import pysam
import primer3

def get_ref_dict(reffilename):
    reffile = pysam.FastaFile(reffilename)
    return {r : reffile.fetch(region=r) for r in reffile.references}

def make_primers_in_vcf(vcffile, referencedict):
    print('\t'.join(['chr', 'pos', 'left', 'right', 'leftpos', 'rightpos']))
    vcfreader = vcf.Reader(filename = vcffile)
    for record in vcfreader:
        res = make_primer_for_record(record, referencedict)
        res = short_result(res)
        if not res['left'] or not res['right']:
            continue
        line = array_result(res)
        print('\t'.join(line))

def make_primer_for_record(record, referencedict):
    chrom = record.CHROM
    pos = record.POS
    seqid = str(chrom) + ":" + str(pos)
    seq_args = {
        'SEQUENCE_ID' : seqid,
        'SEQUENCE_TEMPLATE' : referencedict[chrom],
        'SEQUENCE_TARGET' : [pos, 1]
    }

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

def main():
    refdict = get_ref_dict('/storage/2017-03-15-eucalyptus/e_mel_3/e_mel_3.fa')
    make_primers_in_vcf('/storage/2017-03-15-eucalyptus/reeds_calls/filtered_calls/adams_calls.vcf.gz',refdict)

if __name__ == '__main__':
    main()
