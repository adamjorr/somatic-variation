#!/usr/bin/env python2
import pysam

def load_muts():
	sample_files = ["mut_files/mut_M" + str(k) + ".txt" for k in range(1,9)] #get every file
	muts = []
	for fn in sample_files:
		with open(fn) as fh:
			d = [l.split()[0:1] for l in fh] #a list of every [scaffold,site] in the sample
			muts.append(d)
	return muts #a list of lists of lists

def load_vcf():


def load_all_muts():
	mutfile = "mut_files/mutfile.txt"
	with open(mutfile) as fh:
		for l in fh:
			outline = generate_table_line(l)

def position_depth(samfile, chr, position):
	regionstr = str(chr) + str(position)
	for col in samfile.pileup(regionstr):
		return col.n #just first position in pileup

"""
scaffold, site, original_genotype, mutated_genotype, depth, branch_mutated, samples_mutated, mutation_recovered
"""
def generate_table_line(line, muts, vcf, sam):
	l = line.rstrip().split()
	loc = l[0:1]
	gt = l[2:3]
	togt = l[4:5]
	depth = position_depth(sam, l[0], l[1])
	mutated_samples = ["M" + str(i) for i in range(1,9) if loc in muts[i]]
	

	recovered = False
	if loc in vcfloc:
		if set(vcfgt) == set(togt):
			recovered = True



	return loc[0], loc[1], ''.join(gt), ''.join(togt), ','.join(mutated_samples), 

def main():
    #open all mut files and load into list
    #iterate through full list and generate table

if __name__ == '__main__':
    main()



