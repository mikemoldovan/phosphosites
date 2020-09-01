"""
nonphos_sitelist

Make a list of non-phosphorilated aminoacids in ordered/disordered regions

Input:
1. Proteome FASTA
2. File with disordered regions
3. File with phosphosite list
4. Aminoacid list (default "STY")
5. Outfile name
6. Random var: which percentage of sites should be printed?

Output:
List of non-phosphorilated aminoacids in disordered regions
"""

from optparse import OptionParser
from random import random
from disord_reg_adlib import *


def main(psite_filename,
		 dr_filename,
		 proteome_fasta,
		 aminoacids,
		 outfile_name,
		 string_p):
	
	fastadict = make_fastadict(proteome_fasta)
	dr_dict = build_dr_dict(dr_filename)
	psite_dict = make_psite_dict(psite_filename)

	sitenum_dict = dict()
	psitenum_dict = dict()

	with open(outfile_name, 'w') as outhandle:
		for seq_id in fastadict.keys():
			for dr_crd_pair in dr_dict[seq_id]:
				for i in range(dr_crd_pair[0], dr_crd_pair[1]):
					c = fastadict[seq_id][i]
					if c not in aminoacids:
						continue
					if sitenum_dict.get(c):
						sitenum_dict[c] += 1
					else:
						sitenum_dict[c] = 1 
					if psite_dict.get(seq_id) and psite_dict[seq_id].get(i):
						if psitenum_dict.get(c):
							psitenum_dict[c] += 1
						else:
							psitenum_dict[c] = 1
						continue
					if random() < string_p:
						outhandle.write("{}\t*\t{}\t{}\t*\n".format(seq_id, i, c))
	print sitenum_dict
	print psitenum_dict


parser = OptionParser()
parser.add_option("-p", "--psite_filename", help="File containing a list of phosphosites")
parser.add_option("-d", "--dr_filename", help="File containing disordered (ordered) regions")
parser.add_option("-f", "--proteome_fasta", help="Proteome file")
parser.add_option("-o", "--outfile_name", help="Name of the output file")
parser.add_option("-b", "--string_p", help="Random var: which percentage of sites should be printed? (0-1)")
parser.add_option("-a", "--aminoacids", help="Aminoacids for which a file should be calculated (def = STY)", default="STY")
opt, args = parser.parse_args()


main(opt.psite_filename,
	 opt.dr_filename,
	 opt.proteome_fasta,
	 opt.aminoacids,
	 opt.outfile_name,
	 eval(opt.string_p))

