"""
oma_oggs_spec_filter
Leave in alignments only sequences from selected species

Input:
1. Species list
2. Name of the directory with FASTA sequences
3. Name of the output directory

Output:
Directory with filtered files
"""


import os
from optparse import OptionParser
from Bio import SeqIO


def make_species_dict(species_list_file):
	species_dict = dict()
	with open(species_list_file) as sp_list:
		for s in sp_list:
			s = s.strip()
			s = s.split('_')
			s = ' '.join(s)
			species_dict[s] = True
	return species_dict


def filter_seqs(fafile, fasta_dir, outdir_name, species_dict):
	with open(fasta_dir + fafile) as infile, open(outdir_name + fafile, "w") as outfile:
		fasta_dict = dict()
		for record in SeqIO.parse(infile, "fasta"):
			spec = record.description.split('|')[3].strip().strip('[').strip(']')
			if species_dict.get(spec):
				if fasta_dict.get(spec):
					fasta_dict[spec].append(record)
				else:
					fasta_dict[spec] = [record]

		count = 0
		for spec in fasta_dict.keys():
			for record in fasta_dict[spec]:
				SeqIO.write(record, outfile, "fasta")
			count += 1

	if count == 0:
		os.remove(outdir_name + fafile)


def main(fasta_dir, outdir_name, species_list_file):
	try:
		os.mkdir(outdir_name)
	except:
		pass

	species_dict = make_species_dict(species_list_file)

	for fafile in os.listdir(fasta_dir):
		filter_seqs(fafile, fasta_dir, outdir_name, species_dict)


parser = OptionParser()
parser.add_option("-s", "--species_list_file", help="Species list")
parser.add_option("-i", "--fasta_dir", help="Name of the input directory")
parser.add_option("-o", "--outdir_name", help="Name of the output directory")
opt, args = parser.parse_args()

main(opt.fasta_dir, opt.outdir_name, opt.species_list_file)