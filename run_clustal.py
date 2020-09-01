"""
run_clustal

generate a set of scripts that run clustal on a set of fasta files

Input:
1. Number of scripts
2. Input directory
3. Output directory
4. Script prefix

Output:
Set of scripts
"""


from os import listdir
from optparse import OptionParser


def main(scr_num, fasta_dir, outdir_name, scr_prefix):
	scr_list = []
	for i in range(scr_num):
		scr_file = open(scr_prefix + '_' + str(i) + '.sh', "w")
		scr_list.append(scr_file)
		scr_file.write("#!/bin/bash\n")

	count = 0
	for fasta in listdir(fasta_dir):
		command = "~/tools/clustalo-1.2.4-Ubuntu-x86_64 --in {} --out {}\n".format(fasta_dir + fasta, outdir_name + fasta)
		scr_list[count%scr_num].write(command)
		count += 1

	for scr in scr_list:
		scr.close()


parser = OptionParser()
parser.add_option("-n", "--scr_num", help="Number of scripts")
parser.add_option("-f", "--fasta_dir", help="Directory with FASTA alignments")
parser.add_option("-o", "--outdir_name", help="Output directory name")
parser.add_option("-p", "--scr_prefix", help="Script prefix")
opt, args = parser.parse_args()


main(eval(opt.scr_num), opt.fasta_dir, opt.outdir_name, opt.scr_prefix)