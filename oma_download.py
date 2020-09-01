"""
oma_download
Download OGGs for a SwissProt proteome in FASTA format

Input:
1. File with proteome in FASTA format
2. Name of the output directory
"""

import os
from optparse import OptionParser
from Bio import SeqIO


def main(fasta_proteome, outdir_name):
	try:
		os.mkdir(outdir_name)
	except:
		pass

	os.system("cd {}".format(outdir_name))

	with open(fasta_proteome) as proteome:
		for record in SeqIO.parse(proteome, "fasta"):
			id_str = record.id.split('|')[2]

			os.system('curl "https://omabrowser.org/oma/vps/{}/fasta/" > {}/{}_ogg.fa'.format(id_str, outdir_name, id_str))


parser = OptionParser()
parser.add_option("-i", "--fasta_proteome", help="File with proteome in FASTA format")
parser.add_option("-o", "--outdir_name", help="Name of the output directory")
opt, args = parser.parse_args()

main(opt.fasta_proteome, opt.outdir_name)