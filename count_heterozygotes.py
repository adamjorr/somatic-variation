#!/usr/bin/env python
#python count_heterozygotes.py filein filetype
"""
Deprecated
"""

import sys
import diploidify
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
	'.' : '..',
	'-' : '--'
}

het = ['R','Y','S','W','K','M']

def main():
	filein = sys.argv[1] or sys.stdin
	filetype = sys.argv[2]
	seqs = AlignIO.read(filein,filetype)
	newalignment= [''] * seqs.get_alignment_length()
	for i in range(0,seqs.get_alignment_length()):
		baselist = list(seqs[:,i])
		if not (diploidify.is_valid_site(baselist) and diploidify.is_biallelic(baselist)): continue
		c, v = diploidify.diploidify(baselist)
		if not diploidify.is_single_mutation(c,v): continue
		newalignment = [newalignment[i] + baselist[i] for i in range(len(baselist))]
	newseqobjs = [SeqRecord(Seq(newalignment[i], IUPAC.unambiguous_dna), id=seqs[i].id, description='') for i in range(0,len(seqs))]
	newalnobj = MultipleSeqAlignment(newseqobjs)

	for record in newalnobj:
		print record.id + "\t" + str(len([b for b in record.seq if b in het]))

if __name__ == '__main__':
	main()