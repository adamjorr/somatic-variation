#!usr/bin/env python2

import sys
import argparse
import pysam
import operator

def parse_my_args():
	parser = argparse.ArgumentParser(description="Compare two BAM files built with the same reads and determine where the reads moved to and from.")
	parser.add_argument("-b","--binsize", required=True, help="Size to bin contigs, in Megabases.", type=int)
	parser.add_argument("-n","--maxcontigs", help="The maximum number of contigs to consider; larger contigs are considered first", type=int, default=10)
	parser.add_argument("-o","--out", required=True, help="output file name")
	parser.add_argument("bam1", required=True, help="A BAM file to compare with")
	parser.add_argument("bam2", required=True, help="A BAM file that will be used as a baseline")


	args = parser.parse_args()
	binsize = args.binsize
	maxcontigs = args.maxcontigs
	out = args.out
	otherbam = args.bam1
	basebam = args.bam2

	return (binsize,maxcontigs,out,otherbam,basebam)

def read_bam(file):
	bamob = pysam.AlignmentFile(file,rb)
	bamob.check_index()
	return bamob

def get_chr_tuples(bamob,maxcontigs):
	chrs = bamob.references()
	chr_lengths = bamob.lengths()
	sorted_chrl = sorted(zip(chrs,chr_lengths),key=operator.itemgetter(1))
	return sorted_chrl[:maxcontigs]

def calculate_bins(chr_tuples,binsize):
	bins = []
	for x in chr_tuples:
		start_ends = range(0,x[1],binsize)
		bins.extend([(x[0],start_ends[i],start_ends[i+1]) for i in range(len(start_ends-1))])
		if start_ends[-1] != x[1]:
			bins.extend([(x[0],start_ends[-1],x[1])])

def get_reads_in_bin(bamob,chrm,start,stop):
	bamob.fetch(chrm, start, stop)

def get_unique_reads_in_bin():

def read_movement_value():

def bin_total_movement_value():

def plot_bin_totals():


def main():
	(binsize,maxcontigs,out,otherbam,basebam) = parse_my_args()

if __name__ == '__main__':
	main()
