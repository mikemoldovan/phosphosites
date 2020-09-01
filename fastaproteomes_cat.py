"""
fastaproteomes_cat
Concatenates proteomes in a non-redundant file

Input:
1. List of FASTA files organized by priority
2. Output file name

Output:
Concatenated proteome
"""

from optparse import OptionParser
from Bio import SeqIO

def append_to_fastadict(fafile, fastadict):
	with open(fafile) as infasta:
		for record in SeqIO.parse(infasta, "fasta"):
			rec_id = record.id.split('|')[2]
			if not fastadict.get(rec_id):
				fastadict[rec_id] = record


def main(fastalist, outfasta):
	fastadict = dict()
	for fafile in fastalist.split(','):
		append_to_fastadict(fafile, fastadict)
	with open(outfasta, 'w') as outhandle:
		for k in fastadict.keys():
			SeqIO.write(fastadict[k], outhandle, "fasta")


parser = OptionParser()
parser.add_option("-i", "--fastalist", help="List of FASTA files organized by priority")
parser.add_option("-o", "--outfasta", help="Output file name")
opt, args = parser.parse_args()

main(opt.fastalist, opt.outfasta)