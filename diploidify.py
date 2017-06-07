#!/usr/bin/env python2

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


# There might be a case A/G -> A/T -> A/A = RWA at one site.
# But this is most likely an artifact so we filter it.
def is_valid_site(baselist):
	"""Returns TRUE if we have an iupac entry for all the sites; false otherwise"""
	if False in [b in iupac.keys() for b in baselist]: return False
	return True

def is_biallelic(baselist):
	"""Returns TRUE if list of bases is biallelic; false otherwise"""
	return unique_list_size(baselist) < 3

def unique_list_size(alist):
	"""Returns the size of a unique version of the input list"""
	return len(set(alist))

def non_gap_bases(baselist):
	"""Returns all bases in a list that are not gaps"""
	return [b for b in baselist if b not in gaps]

def is_single_mutation(constantlist, variedlist):
	"""Returns TRUE if only a single mutation occured at this site; false otherwise"""
	cbases = non_gap_bases(constantlist)
	vbases = non_gap_bases(variedlist)
	if unique_list_size(cbases) == 1 and unique_list_size(vbases) == 2: return True
	if unique_list_size(vbases) == 1 and unique_list_size(cbases) == 2: return True 
	if unique_list_size(cbases) == unique_list_size(vbases) == 1: return True
	return False

def diploidify(baselist):
	"""Makes every site 2 sites based on IUPAC codes."""
	constant = [iupac[baselist[0]][0]]
	varied = [iupac[baselist[0]][1]]
	for b in range(1,len(baselist)):
		base = baselist[b]
		if base == '.' or base =='-':
			constant.extend('-')
			varied.extend('-')
		else:
			#We need to be careful about this to account for something like M -> S mutant (A/C -> C/G)
			#So we make a nonvariable haplotype + a variable haplotype
			if constant[0] == iupac[base][0]:
				constant.extend(iupac[base][0])
				varied.extend(iupac[base][1])
			else:
				constant.extend(iupac[base][1])
				varied.extend(iupac[base][0])
	if unique_list_size(constant) != 1:
		dummy = varied
		varied = constant
		constant = dummy
	return constant,varied

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

def remove_nonvariable_sites(baselist):
	"""Takes a baselist and returns an empty list if it is nonvariable, otherwise returns it"""
	if unique_list_size(non_gap_bases(baselist)) == 1:
		return [''] * len(baselist)
	else:
		return baselist

def main():
	parser = argparse.ArgumentParser(description="Create explicit alignment from a diploid sequence that utilizes IUPAC codes")
	parser.add_argument("-i","--input", help="input file (default: stdin)")
	parser.add_argument("-t","--type", help="filetype (default: fasta)", default="fasta")
	parser.add_argument("-o","--output",help="output file (default: stdout)")
	parser.add_argument("-p","--outtype",help="output filetype (default: fasta)", default="fasta")
	parser.add_argument("-v","--variable",help="only output variable sites",action="store_true")
	parser.add_argument("-s","--skip",help="keep single-letter notation but keep other checks",action="store_true")
	args = parser.parse_args()
	filetype = args.type
	filein = args.input or sys.stdin
	fileout = args.output or sys.stdout
	outtype = args.outtype
	variablesites = args.variable
	skip = args.skip

	seqs = AlignIO.read(filein,filetype)
	newalignment= [''] * len(seqs)
	for i in range(0,seqs.get_alignment_length()):
		baselist = list(seqs[:,i])
		if not is_biallelic(baselist): continue
		# if not is_valid_site(baselist): continue
		#Since we have a maximum of 1 mutation, at ambiguous sites we have a constant site and a variable site
		c, v = diploidify(baselist) #turn baselist into 2 baselists, expanding using IUPAC notation
		if not is_single_mutation(c,v): continue
		
		if variablesites:
			c = remove_nonvariable_sites(c)
			v = remove_nonvariable_sites(v)
			if c == v == [''] * len(c): continue
				
		combined = list()
		if skip:
			combined = baselist
		else:
			combined = [c[k] + v[k] for k in range(0,len(c))]
		newalignment = [newalignment[j] + combined[j] for j in range(0,len(combined))]

	newseqobjs = [SeqRecord(Seq(newalignment[l], IUPAC.unambiguous_dna), id=seqs[l].id, description='') for l in range(0,len(seqs))]
	newalnobj = MultipleSeqAlignment(newseqobjs)
	newalnobj = remove_duplicate_seqs(newalnobj)
	AlignIO.write(newalnobj,fileout,outtype)

if __name__ == '__main__':
	main()
