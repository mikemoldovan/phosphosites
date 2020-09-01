"""
godon_run

Prepare scripts for multiple simultaneous godon runs

Input:
1. Directory with alignments
2. Directory with trees
3. Script number
4. Script identifier
5. Name of the output directory

Output:
Set of bash scripts
"""

from optparse import OptionParser
from os import listdir, mkdir

def main(nucl_aldir, tree_dir, scr_num, scr_id, outdir_name):
	try:
		mkdir(outdir_name)
	except:
		pass
	script_files = []
	for i in range(scr_num):
		outfile = open(scr_id + str(i) + '.sh', 'w')
		outfile.write("#!/bin/bash")
		script_files.append(outfile)
	count = 0
	for al in listdir(nucl_aldir):
		alfile_name = nucl_aldir + al
		treefile_name = tree_dir + al.split('.')[0] + ".nwk"
		outfile_name = outdir_name + al.split('.')[0] + ".txt"
#		outfile_name2 = outdir_name + al.split('.')[0] + "2.txt"

		cmd = "~/tools/godon test M8 --codon-omega -p 1 --m0-tree {} {} -o {}\n".format(alfile_name, treefile_name, outfile_name)

		script_files[count%scr_num].write(cmd)
		count += 1
	for scr in script_files:
		scr.close()

#~/tools/godon test M8 --codon-omega --m0-tree human_oggs_no_parals_mammalia_nucl_als2/PCLO_HUMAN_ogg.fa  _PCLO_HUMAN_ogg.nwk

parser = OptionParser()
parser.add_option("-n", "--nucl_al_dir", help="Nucleotide sequence alignment directory")
parser.add_option("-t", "--tree_dir", help="Directory with trees")
parser.add_option("-s", "--scr_num", help="Script number")
parser.add_option("-i", "--scr_id", help="Script identifier")
parser.add_option("-o", "--outdir_name", help="Name of the output directory")
opt, args = parser.parse_args()

main(opt.nucl_al_dir, opt.tree_dir, eval(opt.scr_num), opt.scr_id, opt.outdir_name)