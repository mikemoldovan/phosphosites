"""
prepare_als_for_paml

Convert alignments in a given directory to PHYLIP and reduce FASTA description
to species names

Input:
1. Directory with FASTA alignments
2. Output directory name

Output:
Directory with processed alignments
"""


from os import listdir, mkdir, remove
from optparse import OptionParser
from Bio import SeqIO, AlignIO


def build_spec_code_dict(spec_code_table):
	spec_code_dict = dict()
	with open(spec_code_table) as sct:
		for s in sct:
			s = s.strip().split()
			spec_code_dict[s[0]] = s[1]
	return spec_code_dict


def print_phylip(infasta, outfile):
	seq_dict = dict()
	length = 0
	seq_num = 0
	for record in SeqIO.parse(infasta, "fasta"):
		seq_num += 1
		seq = ""
		length = 0
		for c in record.seq:
			seq += c
			length += 1
		seq_dict[record.id] = seq
	outfile.write(" {}  {}\n\n".format(str(seq_num), str(length)))
	for k in seq_dict.keys():
		outfile.write(k + ' '*(11 - len(k)))
		count = 0
		for c in seq_dict[k]:
			count += 1
			outfile.write(c)
			if count % 10 == 0 and count != length:
				outfile.write(' ')
		outfile.write('\n')


def process_al(fafile_name, fasta_dir, outdir_name, spec_code_dict):
	phylname = '.'.join(fafile_name.split('.')[:-1]) + ".phy"
	
	with open(fasta_dir + fafile_name) as infasta, open("temp.fa", "w") as temp:
		for record in SeqIO.parse(infasta, "fasta"):
			spec = record.description.split('|')[3].strip().strip('[').strip(']')
			spec = '_'.join(spec.split())
			seq = ""
			for c in record.seq:
				seq += c
			temp.write('>' + spec_code_dict[spec] + '\n')
			temp.write(seq + '\n')

	with open("temp.fa") as temp, open(outdir_name + phylname, "w") as phylip_al:
		print_phylip(temp, phylip_al)
#		als = AlignIO.parse(temp, "fasta")
#		AlignIO.write(als, phylip_al, "phylip")
#		phylip_al.write(' ')

#	remove("temp.fa")


def main(fasta_dir, outdir_name, spec_code_table):
	spec_code_dict = build_spec_code_dict(spec_code_table)

	try:
		mkdir(outdir_name)
	except:
		pass

	for fafile_name in listdir(fasta_dir):
		process_al(fafile_name, fasta_dir, outdir_name, spec_code_dict)
		


parser = OptionParser()
parser.add_option("-f", "--fasta_dir", help="Directory with FASTA alignments")
parser.add_option("-o", "--outdir_name", help="Output directory name")
parser.add_option("-s", "--spec_code_table", help="Table with species codes (prepare_tree_for_paml.py output)")
opt, args = parser.parse_args()


main(opt.fasta_dir, opt.outdir_name, opt.spec_code_table)