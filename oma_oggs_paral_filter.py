"""
oma_oggs_paral_filter
Filter out all paralogs from fasta files

Input:
1. Name of the directory with FASTA sequences
2. Name of the output directory

Output:
Directory with filtered files
"""

import os
from optparse import OptionParser
from Bio import SeqIO


def filter_seqs(fafile, fasta_dir, outdir_name):
	with open(fasta_dir + fafile) as infile, open(outdir_name + fafile, "w") as outfile:
		spec_dict = dict()
		for record in SeqIO.parse(infile, "fasta"):
			spec = record.description.split('|')[3].strip()
			if spec_dict.get(spec):
				spec_dict[spec].append(record)
			else:
				spec_dict[spec] = [record]

		count = 0
		for spec in spec_dict.keys():
			if len(spec_dict[spec]) > 1:
				continue
			SeqIO.write(spec_dict[spec][0], outfile, "fasta")
			count += 1

	if count == 0:
		os.remove(outdir_name + fafile)


def main(fasta_dir, outdir_name):
	try:
		os.mkdir(outdir_name)
	except:
		pass

	for fafile in os.listdir(fasta_dir):
		filter_seqs(fafile, fasta_dir, outdir_name)


parser = OptionParser()
parser.add_option("-i", "--fasta_dir", help="Name of the input directory")
parser.add_option("-o", "--outdir_name", help="Name of the output directory")
opt, args = parser.parse_args()

main(opt.fasta_dir, opt.outdir_name)

