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
	dipl = [iupac[b.upper()] for b in baselist]
	return create_cv(dipl)

def create_cv(bases):
	"""Turns a list of diploid bases into 2 lists, attempting to create one constant and one varied (but no promises)"""
	constant = [bases[0][0]]
	varied = [bases[0][1]]
	for b in bases[1:]:
		b.replace('.','-')# standardize gaps
		if b[0] == constant[0]:
			constant.extend(b[0])
			varied.extend(b[1])
		else:
			constant.extend(b[1])
			varied.extend(b[0])
	if unique_list_size(constant) != 1:
		constant, varied = varied, constant #switch variables
	return constant, varied

def filter_cv(c, v, args):
	"""Filter constant and variable sites according to arguments"""
	if not is_single_mutation(c,v): return None, None
	if args.variable:
		c = remove_nonvariable_sites(c)
		v = remove_nonvariable_sites(v)
		if c == v == None: return None, None
	return c, v


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
	"""Takes a baselist and returns None if it is nonvariable, otherwise returns it"""
	if unique_list_size(non_gap_bases(baselist)) == 1:
		return None
	else:
		return baselist

def combine_cv(c,v):
	if c and v:
		return [c[k] + v[k] for k in range(0,len(c))]
	else:
		return (v if v else c)

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
		if not is_biallelic(baselist): continue
		#Since we have a maximum of 1 mutation, at ambiguous sites we have a constant site and a variable site
		c, v = diploidify(baselist) #turn baselist into 2 baselists, expanding using IUPAC notation
		# if not is_single_mutation(c,v): continue
		
		# if variablesites:
		# 	c = remove_nonvariable_sites(c)
		# 	v = remove_nonvariable_sites(v)
		# 	if c == v == [''] * len(c): continue
		c, v = filter_cv(c,v,args)
		if c == v == None: continue		
		combined = list()
		if skip:
			combined = baselist
		else:
			combined = combine_cv(c,v)
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
				gts.append(''.join(bases.replace(sepchar,''))) #sort to ensure A/T is not different from T/A
			else:
				gts.append(None)
				break
		if gts[-1] is None: continue
		c,v = create_cv(gts)
		c,v = filter_cv(c,v,args)
		if c == v == None: continue
		if ot == "vcf":
			vcf_writer.write_record(record)
		else:
			combined = combine_cv(c,v)
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


def argparser():
	parser = argparse.ArgumentParser(description=
		"""
		Create explicit alignment from a diploid sequence that utilizes IUPAC codes.
		Will remove sequences that are identical to previous sequences in the alignment.
		Will also remove any site that is not biallelic.
		Use "-t vcf" to perform these checks and remove nonvariable sites from a vcf.
		In VCF mode, this program will only output a vcf.
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

def main():
	args = argparser()
	if args.type != "vcf":
		filter_alignment(args)
	else:
		filter_vcf(args)

if __name__ == '__main__':
	main()
