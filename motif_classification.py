"""
motif_classification

Employs a phosphosite classification approach proposed in:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1785252/

Input:
1. Phosphosite table
2. Proteome file

Output:
Phosphosite table with marked phosphosites
"""

from optparse import OptionParser
from disord_reg_adlib import *

def is_Y(seq, i):
	if seq[i] == 'Y':
		return 'Y'

def is_P(seq, i):
	try:
		if seq[i+1] == 'P':
			return 'P'
	except IndexError:
		return None

def is_A1(seq, i):
	try:
		count = 0
		for j in range(1, 7):
			if seq[i + j] in ('E', 'D'):
				count += 1
		if count >= 5:
			return 'A'
	except IndexError:
		return None

def is_B1(seq, i):
	if i < 3:
		return None
	try:
		if seq[i - 3] in ('R', 'K'):
			return 'B'
	except IndexError:
		return None

def is_A2(seq, i):
	try:
		for j in range(1, 4):
			if seq[i + j] in ('E', 'D'):
				return 'A'
	except IndexError:
		return None

def is_B2(seq, i):
	try:
		count = 0
		for j in range(1, 7):
			if i < j:
				break
			if seq[i - j] in ('R', 'K'):
				count += 1
		if count >= 2:
			return 'B'
	except IndexError:
		return None

def is_other(seq, i):
	return 'O'


def main(psite_filename, proteome_fasta):
#	print psite_filename, proteome_fasta
	psite_dict = make_psite_dict2(psite_filename)
	proteome_fastadict = make_fastadict(proteome_fasta)
	testlist = [is_Y, is_P, is_A1, is_B1, is_A2, is_B2, is_other]

	for uniprot_id in psite_dict.keys():
		if not proteome_fastadict.get(uniprot_id):
			continue
		seq = proteome_fastadict[uniprot_id]
		for crd in psite_dict[uniprot_id].keys():
			for test in testlist:
				testres = test(seq, crd)
				if testres != None:
					outline = psite_dict[uniprot_id][crd]
					outline.append(testres)
					outline = '\t'.join(outline)
					print outline
					break


parser = OptionParser()
parser.add_option("-p", "--psite_filename", help="Phosphosite table")
parser.add_option("-d", "--proteome_fasta", help="Proteome file (optonal, if there is no table with site numbers)", default = "_")
opt, args = parser.parse_args()


main(opt.psite_filename, opt.proteome_fasta)
