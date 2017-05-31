#!/usr/bin/env python3

from __future__ import print_function
import argparse
import screed
import sys
import khmer


def output_single(read):
    if hasattr(read, 'quality'):
        return "@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.quality)
    else:
        return ">%s\n%s\n" % (read.name, read.sequence)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--min-coverage', type=int, default=None)
    parser.add_argument('-M', '--max-coverage', type=int, default=None)
    parser.add_argument('input_count_graph')
    parser.add_argument('input_readfile')
    parser.add_argument('input_pairfile')
    parser.add_argument('output_readfile')
    parser.add_argument('output_pairfile')
    parser.add_argument('output_singlefile', default=None, nargs='?')
    args = parser.parse_args()

    print('min_coverage: %s' % args.min_coverage, file=sys.stderr)
    print('max_coverage: %s' % args.max_coverage, file=sys.stderr)

    if not (args.min_coverage or args.max_coverage):
        print("neither min nor max coverage specified!? exiting!", file=sys.stderr)
        sys.exit(1)

    if args.min_coverage and args.max_coverage and \
       args.max_coverage < args.min_coverage:
        print("min_coverage > max_coverage!? exiting!", file=sys.stderr)
        sys.exit(1)

    htable = khmer.load_countgraph(args.input_count_graph)
    output_file = args.output_readfile
    output_fp = open(output_file, 'w')
    output_pairfile = args.output_pairfile
    output_pfp = open(output_pairfile,'w')

    if args.output_singlefile:
        output_singlefile = args.output_singlefile
        output_sfp = open(output_singlefile,'w')

    n_kept = 0
    n = 0
    pair_iter = iter(screed.open(args.input_pairfile))

    for n, record in enumerate(screed.open(args.input_readfile)):
        if n % 100000 == 0:
            print('...', n, n_kept, file=sys.stderr)

        seq = record.sequence.upper()
        seq = seq.replace('N', 'A')
        pair_record = next(pair_iter)
        pseq = pair_record.sequence.upper()
        pseq = pseq.replace('N','A')

        try:
            med, _, _ = htable.get_median_count(seq)
            pmed, _, _ = htable.get_median_count(pseq)
        except ValueError:
            continue

        keep = True
        pkeep = True
        if args.min_coverage:
            if med < args.min_coverage:
                keep = False
            if pmed < args.min_coverage:
                pkeep = False

        if args.max_coverage:
            if med > args.max_coverage:
                keep = False
            if pmed > args.max_coverage:
                pkeep = False

        if keep and pkeep:
            n_kept += 2
            output_fp.write(output_single(record))
            output_pfp.write(output_single(pair_record))

        elif args.output_singlefile:
            if keep:
                n_kept += 1
                output_sfp.write(output_single(record))
            elif pkeep:
                n_kept += 1
                output_sfp.write(output_single(pair_record))

    print('consumed %d reads; kept %d' % (2 * (n + 1), n_kept), file=sys.stderr)

if __name__ == '__main__':
    main()
