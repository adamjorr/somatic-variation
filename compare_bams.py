#!usr/bin/env python2
#WIP to visualize movement of reads between two BAM files created with the same set of reads.

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
	'''read a file and return a pysam bam object'''
	bamob = pysam.AlignmentFile(file,rb)
	bamob.check_index()
	return bamob

def get_chr_tuples(bamob,maxcontigs):
	'''returns a list of tuples, where each tuple is (chromosome, length of chromosome)'''
	chrs = bamob.references()
	chr_lengths = bamob.lengths()
	sorted_chrl = sorted(zip(chrs,chr_lengths),key=operator.itemgetter(1))
	return sorted_chrl[:maxcontigs]

def calculate_bins(chr_tuples,binsize):
	'''returns a list of tuples, where each tuple is (chromosome, start, stop)'''
	bins = []
	for x in chr_tuples:
		start_ends = range(0,x[1],binsize)
		bins.extend([(x[0],start_ends[i],start_ends[i+1]) for i in range(len(start_ends-1))])
		if start_ends[-1] != x[1]:
			bins.extend([(x[0],start_ends[-1],x[1])])
	return bins

def get_reads_in_bin(bamob,chrm,start,stop):
	'''returns dictionary of read names in region'''
	d = dict()
	for read in bamob.fetch(chrm, start, stop):
		#get a read (pysam.AlignedSegment) in region and...
		d[read.query_name] = 0

def reads_added_to_bin(first, second):
	'''takes dictionaries of reads and returns a list of reads present in second but not first'''
	[k for k in second.keys() if k not in first]

def reads_removed_from_bin(first, second):
	'''takes dictionaries of reads and returns a list of reads present in first but not second'''
	[k for k in first.keys() if k not in second]

def bin_movement_values(bamob,maxcontigs,binsize):
	d = dict()
	chrs = get_chr_tuples(bamob,maxcontigs)
	bins = calculate_bins(chrs,binsize)
	for b in bins:
		d[b] = get_reads_in_bin(bamob,b[0],b[1],b[2])

def compare_bams(bam1,bam2,maxcontigs,binsize):
	bam_obj1 = read_bam(bam1)
	bam_obj2 = read_bam(bam2)
	dictionary1 = bin_movement_values(bam_obj1,maxcontigs,binsize)
	dictionary2 = bin_movement_values(bam_obj2,maxcontigs,binsize)
	return (dictionary1, dictionary2)

def print_bin_totals(dict1,dict2,out):
	bins = list()
	added = list()
	removed = list()
	for k in dict1.keys():
		b.append(k)
		added.append(reads_added_to_bin(dict1[k],dict2[k]))
		removed.append(reads_removed_from_bin(dict1[k],dict2[k]))

	with open(out,'wb') as f:
		f.write("\t".join("bin","added","removed"))
		for b in bins:
			f.write("\t".join(b,len(added),len(removed)))

def main():
	(binsize,maxcontigs,out,otherbam,basebam) = parse_my_args()
	dict1, dict2 = compare_bams(basebam,otherbam,maxcontigs,binsize)
	print_bin_totals(dict1,dict2,out)

if __name__ == '__main__':
	main()
