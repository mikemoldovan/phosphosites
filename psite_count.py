"""
psite_count

Count phosphosites remeined after filtering procedures

Input:
1. Phosphosite table
2. Directory with alignments/OGG fastas

Output:
number of remained phosphosites
"""

from os import listdir
from optparse import OptionParser
from Bio import SeqIO


def make_psite_dict(psite_table):
	psite_dict = dict()
	with open(psite_table) as psites:
		for s in psites:
			s = s.strip().split()
			if s[0] == "UniProtID":
				continue
			if psite_dict.get(s[0]):
				psite_dict[s[0]] += 1
			else:
				psite_dict[s[0]] = 1
	return psite_dict


def infasta(UniProtID, fafile):
	with open(fafile) as infasta:
		for record in SeqIO.parse(infasta, "fasta"):
			_id = record.description.split('|')[1].strip()
			if _id == UniProtID:
				return True
	return False


def main(psite_table, fasta_dir):
	count = 0
	psite_dict = make_psite_dict(psite_table)

	for fasta in listdir(fasta_dir):
		UniProtID = fasta[:-7]
		if infasta(UniProtID, fasta_dir + fasta):
			if psite_dict.get(UniProtID):
				count += psite_dict[UniProtID]
	print count


parser = OptionParser()
parser.add_option("-p", "--psite_table", help="Phosphosite table")
parser.add_option("-i", "--fasta_dir", help="Directory with alignments/OGG fastas")
opt, args = parser.parse_args()

main(opt.psite_table, opt.fasta_dir)