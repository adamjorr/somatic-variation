#!usr/bin/env python
#python supported_sites.py

"""
Prints all the splits in the tree and the number of SNPs in the alignment that support them.

A human-readable report will be printed to stdout.
A tab-separated table can be printed to a file specified by -o
A tab-separated table with the number of homozygous sites for each taxa can be printed to a file specified by -z
If -o or -z are STDOUT, they will also be printed to stdout. Useful with -n.
"""

from Bio import Phylo
from Bio import AlignIO
from Bio import SeqIO
import string
import sys
import argparse

iupac_het = set(['R','Y','S','W','K','M'])
iupac_homo = set(['A','C','T','G'])



def parse_my_args():
	parser = argparse.ArgumentParser(description="Find the number of sites that support each split in a given tree")
	parser.add_argument("-t","--tree", required=True, help="newick tree file", metavar='TREE.nwk')
	parser.add_argument("-a","--alignment", required=True, help="alignment file")
	parser.add_argument("-f","--filetype", help="alignment filetype (default: phylip-relaxed)", default="phylip-relaxed")
	parser.add_argument("-s","--singletons",help="use if you want the number of singleton sites as well",action='store_true')
	parser.add_argument("-o","--out", help="optional file to print support table to, can be STDOUT")
	parser.add_argument("-z","--homozygotes", help="optional file to print table with number of homozygous sites, can be STDOUT")
	parser.add_argument("-n","--noreport",help="doesn't output the human-readable report to STDOUT.",action='store_false')
	args = parser.parse_args()
	alignin = args.alignment
	filetype = args.filetype
	treein = args.tree
	singletons = args.singletons
	out = args.out
	hz = args.homozygotes
	report = args.noreport
	return (alignin, filetype, treein, singletons,out,hz,report)

def get_splits(tree, singletons):
	#Return taxa that compose each split as a list of tuples.
	taxanames = [t.name for t in tree.get_terminals()]
	nonterminals = tree.get_nonterminals()
	splits = []
	if singletons: [splits.append(([t],[o for o in taxanames if o != t])) for t in taxanames]
	for interior in nonterminals:
		ingroup = [t.name for t in interior.get_terminals()]
		outgroup = [t for t in taxanames if t not in ingroup]
		if 1 < len(ingroup) < len(taxanames):
			splits.append((ingroup,outgroup))
	return splits

def base_set(alignment, i):
	bases = [b for b in alignment[:,i] if b not in ['-','.','N']]
	return set(bases)

def site_split(alignment, i):
	bases = base_set(alignment,i)
	if len(bases) > 2: raise ValueError(
		'This sequence is not biallelic! At site ' + str(i) + '!'
		)
	if len(bases) == 0: raise ValueError(
		'This sequence has a site with only gaps!'
		)
	bases = list(bases)
	taxa1 = [t.id for t in alignment if t.seq[i] == bases[0]]
	taxa2 = [t.id for t in alignment if t.seq[i] == bases[1]] if len(bases) > 1 else []
	return (taxa1, taxa2)

def sites_zygosity(alignment, i):
	bases = base_set(alignment,i)
	bases = list(bases)
	if len(bases) == 1:
		bases.extend([None])
	if bases[0] in iupac_het:
		homozygous_base = bases[1]
		heterozygous_base = bases[0]
	else:
		homozygous_base = bases[0]
		heterozygous_base = bases[1]
	homozygous = [t.id for t in alignment if t.seq[i] == homozygous_base]
	heterozygous = [t.id for t in alignment if t.seq[i] == heterozygous_base]
	return (homozygous, heterozygous)

def num_sites(alignment):
	return 

def pct_homozygous(alignment):
	num_hom = {record.id : len([b for b in record.seq if b in iupac_homo]) for record in alignment}
	all_sts = {record.id : len([b for b in record.seq if b not in ['-','.','N']]) for record in alignment}
	return {k:float(v)/all_sts[k] for k,v in num_hom.iteritems()}

def most_homozygous(alignment):
	"""
	Finds the taxon with the most homozygous sites
	"""
	pct_homo = pct_homozygous(alignment)

	return max(pct_homo, key=lambda x : pct_homo[x])

def split_support(splits,taxa1,taxa2):
	#Return a list of splits this site supports
	supporting = list()
	for split in splits:
		split1 = set(split[0])
		split2 = set(split[1])
		if taxa1.issubset(split1) and taxa2.issubset(split2):
			supporting.append(split)
		elif taxa1.issubset(split2) and taxa2.issubset(split1):
			supporting.append(split)
	return supporting

def zygosity_support(splits,homozygous,heterozygous):

def antizygosity_support(splits,homozygous,heterozygous):

def incompatible_split(splits,homozygous,heterozygous):

def match_alignment(splits, alignment):
	root_taxon = most_homozygous(alignment)
	counters = {str(split):[0]*4 for split in splits}
	for i in range(alignment.get_alignment_length()):
		(taxa1, taxa2) = site_split(alignment, i)
		(homozygous, heterozygous) = sites_zygosity(alignment,i)
		taxa1 = set(taxa1)
		taxa2 = set(taxa2)
		homozygous = set(homozygous)
		heterozygous = set(heterozygous)

		supporting = split_support(splits,taxa1,taxa2)
		for split in supporting: counters[str(split)][0] += float(1)/len(supporting)

		#Possibly break this into its own function
		for split in splits:
			split1 = set(split[0])
			split2 = set(split[1])
			left_split = split1 if root_taxon in split1 else split2
			right_split = split1 if root_taxon not in split1 else split2

			#Zygosity Support
			if homozygous.issubset(left_split) and heterozygous.issubset(right_split):
				counters[str(split)][1] += 1

			#Anti-zygosity Support
			if homozygous.issubset(right_split) and heterozygous.issubset(left_split):
				counters[str(split)][2] += 1

			#Incompatible Splits
			if len(split1 & heterozygous) > 0 and len(split1 & homozygous) > 0:
				if len(split2 & heterozygous) > 0 and len(split2 & homozygous) > 0:
					counters[str(split)][3] += 1

	return counters

def print_report(counters,alnlen):
	print '\nAlignment length: ' + str(alnlen) + "\n"
	for k,v in counters.iteritems():
		w = [float(w)/alnlen for w in v]
		print str(k) + "\t" + str(v) + " / " + "{0:.0%} {1:.0%} {2:.0%} {3:.0%}".format(w[0],w[1],w[2],w[3])
	print ''

def write_support_table(f,counters,alnlen):
	header = '\t'.join(['split','support','zygosity_support','anti_zygosity','incompatible'])
	f.write(header + '\n')
	for k,v in counters.iteritems():
		v = [float(w)/alnlen*100 for w in v]
		tableval = k[1:-1].split('[')[1].split(']')[0].replace(' ','').replace('\'','') +"\t{0:.0f}\t{1:.0f}\t{2:.0f}\t{3:.0f}".format(v[0],v[1],v[2],v[3])
		f.write(tableval + '\n')

def write_hom_table(f,alignment):
	num_homo = pct_homozygous(alignment)
	header = '\t'.join(['taxon','percent_homozygous'])
	f.write(header + '\n')
	for k,v in num_homo.iteritems():
		v = v*100
		f.write(k + '\t{0:.0f}'.format(v) + '\n')

def main():
	(alignin, filetype, treein, singletons,out,hz,report) = parse_my_args()
	alignment = AlignIO.read(alignin,filetype)
	tree = Phylo.read(treein,'newick')
	splits = get_splits(tree,singletons)
	counters = match_alignment(splits, alignment)
	alnlen = alignment.get_alignment_length()

	if report:
		print_report(counters,alnlen)
	
	if out != None:
		f = open(out,'w') if out.lower() != 'stdout' else sys.stdout
		write_support_table(f,counters,alnlen)
		if f != sys.stdout:
			f.close()
		else: f.write('\n')

	if hz != None:
		if hz.lower() != 'stdout':
			g = open(hz,'a') if hz == out else open(hz,'w')
		else:
			g = sys.stdout
		write_hom_table(g,alignment)
		if g != sys.stdout:
			g.close()
		else: g.write('\n')


if __name__ == '__main__':
	main()