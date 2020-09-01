"""
build_and_trim_alignments

runs build_and_trim_alignment.py on a directory

Input:
1. Protein sequence alignment directory
2. Nucleotide sequence alignment directory
3. Name of the output directory
"""

from optparse import OptionParser
from os import system, listdir, mkdir

def main(protseq_dir, nucl_seq_dir, out_dir):
	try:
		mkdir(out_dir)
	except:
		pass
	for filename in listdir(nucl_seq_dir):
		system("python ~/phosphosites/scripts/build_and_trim_alignment.py -n {} -p {} -o {}".format(nucl_seq_dir + filename, protseq_dir + filename, out_dir + filename))

parser = OptionParser()
parser.add_option("-n", "--nucl_seq_dir", help="Nucleotide sequence directory")
parser.add_option("-p", "--protseq_dir", help="Corresponding protein sequence directory")
parser.add_option("-o", "--out_dir", help="Name of the output directory")
opt, args = parser.parse_args()


main(opt.protseq_dir, opt.nucl_seq_dir, opt.out_dir)