"""
nucl_seq_oggs

Infer nucleotide sequences fpr corresponding OMA protein files from a given directory
using the nucleotide sequence file obtained from https://omabrowser.org/oma/current/

Input:
1. Nucleotide sequence file
2. Protein sequence directory
3. Name of the output directory

Output:
Directory containing nucleotide sequences matching those from protein files
"""

from os import listdir, mkdir
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna
from Bio import SeqIO
from optparse import OptionParser

def fastaparse(fafile):
	firstline = True
	with open(fafile) as infasta:
		for s in infasta:
			if s.startswith('>'):
				if not firstline:
					yield seq_id, seq
					seq_id = s.strip().split()[1]
					seq = ""
				else:
					seq_id = s.strip().split()[1]
					firstline = False
					seq = ""
			else:
				seq += s.strip()
	yield seq_id, seq 


def makeprotdict(protfasta_dir):
	protdict = dict()
	ogg_dict = dict()
	for protfasta in listdir(protfasta_dir):
		filename = protfasta_dir + protfasta
		print(filename)
		ogg_dict[protfasta] = dict()
		with open(filename) as inhandle:
			for record in SeqIO.parse(inhandle, "fasta"):
				seq = ""
				for c in record.seq:
					seq += c
				protdict[record.id] = seq
				ogg_dict[protfasta][record.id] = True
#		break
	print(protdict, ogg_dict)
	return protdict, ogg_dict


def main(nucl_seq_file, protfasta_dir, outdir_name):
	try:
		mkdir(outdir_name)
	except:
		pass
	protdict, ogg_dict = makeprotdict(protfasta_dir)
	nucl_dict = dict()
	for seq_id, seq in fastaparse(nucl_seq_file):
		if protdict.get(seq_id):
			nucl_dict[seq_id] = seq
	print(nucl_dict)
	for ogg_id in ogg_dict.keys():
		outhandle = open(outdir_name + ogg_id, 'w')
		for seq_id in ogg_dict[ogg_id].keys():
			outhandle.write('>' + seq_id + '\n')
#			nucl_seq = nucl_dict[seq_id]
#			rna = Seq(nucl_seq, generic_rna)
#			trans_seq = ""
#			for c in rna.translate():
#				trans_seq += c
#			if trans_seq[:-1] == protdict[seq_id]:
			outhandle.write(nucl_dict[seq_id] + '\n')
		outhandle.close()


parser = OptionParser()
parser.add_option("-n", "--nucl_seq_file", help="Nucleotide sequence file")
parser.add_option("-p", "--protfasta_dir", help="Protein sequence directory")
parser.add_option("-o", "--outdir_name", help="Name of the output directory")
opt, args = parser.parse_args()

main(opt.nucl_seq_file, opt.protfasta_dir, opt.outdir_name)


