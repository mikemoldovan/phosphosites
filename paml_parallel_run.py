"""
paml_parallel_run

Run PAML multiple times and via multiple scripts:
Produces a set of .sh scripts for "multithread" PAML run with paml_run.py script

Input:
1. Output directories ID
2. Number of scripts
3. Directory containing PHYLIP sequences for PAML
4. Directory containing NEWICK trees for paml
5. PAML config file (FULL PATH: python)
/mnt/mapr/user/mmoldovan/phosphosites/data/sample_paml_config.ctl

Output:
scripts for the PAML run
"""


from optparse import OptionParser
from os import listdir


def main(outdir_id, scr_num, phylip_dir, nwk_dir, paml_ctl):
	dirlist = []
	for i in range(scr_num):
		dirlist.append(outdir_id + '_' + str(i) + '/')

	scr_list = []
	for i in range(scr_num):
		scr = open("codeml_run_" + str(i) + ".sh", "w")
		scr_list.append(scr)
		scr.write("#!/bin/bash\n#$ -pe smp 1\n#$ -cwd\n")

	count = 0
	for nw_tree_fname in listdir(phylip_dir):
		count += 1
		job_id = nw_tree_fname.split('/')[-1][:-8]
		outdir_name = dirlist[count % scr_num]
		phy_align_fname = phylip_dir + job_id + "_ogg.phy"
		nw_tree_fname = nwk_dir + job_id + "_ogg.nwk"
		scr_list[count % scr_num].write("python ~/phosphosites/scripts/paml_run.py -n {} -a {} -o {} -c {}\n".format(nw_tree_fname,
																												  phy_align_fname,
																												  outdir_name,
																												  paml_ctl))

	for i in range(scr_num):
		scr_list[i].close()


parser = OptionParser()
parser.add_option("-n", "--scr_num", help="Number of scripts")
parser.add_option("-a", "--phylip_dir", help="Directory containing PHYLIP sequences for PAML (FULL PATH: unix)")
parser.add_option("-o", "--outdir_id", help="Output directories ID (FULL PATH: unix)")
parser.add_option("-w", "--nwk_dir", help="Directory containing NEWICK trees for paml (FULL PATH: unix)")
parser.add_option("-c", "--paml_ctl", help="PAML config file (FULL PATH: python)")
opt, args = parser.parse_args()


main(opt.outdir_id, eval(opt.scr_num), opt.phylip_dir, opt.nwk_dir, opt.paml_ctl)
