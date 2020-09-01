"""
paml_run

Takes a file with a tree and, a file with an alignment and a PAML config file

Runs PAML and produces two files: a NEWICK file containing a tree and a FASTA with
both ancestral and contemporary sequences

Input:
1. NEWICK tree file name (FULL PATH: unix)
2. PHYLIP alignment file name (FULL PATH: unix)
3. PAML config file (FULL PATH: python)
/mnt/mapr/user/mmoldovan/phosphosites/data/sample_paml_config.ctl
4. Directory for the output

Output:
1. FASTA alignment file for all sequences
2. PHYLIP tree file containing a tree
"""


from optparse import OptionParser
from os import remove, chdir, system, getcwd, mkdir


def run_codeml(nw_tree_fname, phy_align_fname, paml_config):
	config_str = ""
	with open(paml_config) as p_c:
		for s in p_c:
			config_str += s

	config_str = config_str.replace("[1]", phy_align_fname)
	config_str = config_str.replace("[2]", nw_tree_fname)

	with open("paml_config.ctl", "w") as ctl:
		ctl.write(config_str)

	system("~/tools/paml4.9h/bin/codeml paml_config.ctl")


def codeml_res_parse(job_id):
	inseqblock = False
	count_seqblock = 0
	with open("rst") as rst, open(job_id + ".nwk", "w") as treefile, open(job_id + ".fa", "w") as alfile:
		for s in rst:
			if "tree with node labels for Rod Page's TreeView" in s:
				treefile.write(next(rst))
			if "List of extant and reconstructed sequences" in s:
				inseqblock = True
			if "Overall accuracy" in s:
				inseqblock = False
			if inseqblock:
				s = s.strip().split()
				if len(s) >= 2:
					if s[0].isdigit() and not s[1].isdigit():
						count_seqblock += 1
						alfile.write('>' + str(count_seqblock) + '_' + s[0] + '\n')
						alfile.write(''.join(s[1:]) + '\n')
					elif s[0] == 'node':
						alfile.write('>' + s[1][1:] + '\n')
						alfile.write(''.join(s[2:]) + '\n')
	remove("test.mlc")
	remove("rub")
	remove("rst1")
	remove("rst")
	remove("rates")
	remove("lnf")


def main(nw_tree_fname, phy_align_fname, paml_config, output_dir):
	try:
		mkdir(output_dir)
	except:
		pass

	cwd = getcwd()

	chdir(output_dir)
	job_id = nw_tree_fname.split('/')[-1][:-8]
	run_codeml(nw_tree_fname, phy_align_fname, paml_config)
	codeml_res_parse(job_id)
	chdir(cwd)


parser = OptionParser()
parser.add_option("-n", "--nw_tree_fname", help="NEWICK tree file name (FULL PATH: unix)")
parser.add_option("-a", "--phy_align_fname", help="PHYLIP alignment file name (FULL PATH: unix)")
parser.add_option("-o", "--output_dir", help="Directory for the output")
parser.add_option("-c", "--paml_config", help="PAML config file (FULL PATH: python)")
opt, args = parser.parse_args()


main(opt.nw_tree_fname, opt.phy_align_fname, opt.paml_config, opt.output_dir)