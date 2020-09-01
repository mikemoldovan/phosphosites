"""
prepare_trees_for_paml

For each alignment file prune the "big tree" leaving only
leaves corresponding to organisms presented in alignment files.

Resulting trees can be used as PAML input

Input:
1. Directory with alignment files
2. Big tree file
3. Output directory name

Output:
Directory with pruned trees
"""


from optparse import OptionParser
from Bio import SeqIO
#from copy import deepcopy
from os import listdir, mkdir, remove, getcwd, chdir

#cwd = getcwd()
#print getcwd()
#chdir("/mnt/mapr/user/mmoldovan/tools/ete_master")
#print getcwd()
#from ete3 import Tree
#chdir(cwd)


def make_leaf_list(phylip_input_name, phylip_input_dir):
	leaf_list = []
	leaf_num = 0
	with open(phylip_input_dir + phylip_input_name) as infile:
		for s in infile:
			s = s.strip().split()
			if len(s) < 2:
				continue
			if s[1].isdigit():
				continue
			leaf_list.append(s[0])
			leaf_num += 1
#	print phylip_input_name, leaf_list
	return leaf_list, leaf_num


def print_pruned(pruned_tree, leaf_num, tree_filename):
	pruned_tree.write(format=1, outfile="temp.nw")
	with open(tree_filename, "w") as tree_outfile, open("temp.nw") as temp:
		tree_outfile.write(" {}  1\n\n".format(str(leaf_num)))
		for s in temp:
			tree_outfile.write(s)
	remove("temp.nw")


def main(alfile_dir, bigtree_filename, outdir_name):
	try:
		mkdir(outdir_name)
	except:
		pass
	bigtree_ref = Tree(bigtree_filename)
	for alfile in listdir(alfile_dir):
		bigtree = deepcopy(bigtree_ref)
		fasta_input_name = alfile_dir + alfile
		leaf_list, leaf_num = make_leaf_list(alfile, alfile_dir)
		bigtree.prune(leaf_list)
		tree_filename = outdir_name + '.'.join(alfile.split('.')[:-1]) + ".nwk"
		print_pruned(bigtree, leaf_num, tree_filename)
		break


parser = OptionParser()
parser.add_option("-a", "--alfile_dir", help="Directory with alignment files")
parser.add_option("-b", "--bigtree_filename", help="Big tree file")
parser.add_option("-o", "--outdir_name", help="Output directory name")
opt, args = parser.parse_args()


main(opt.alfile_dir, opt.bigtree_filename, opt.outdir_name)

#prepare_trees_for_paml.main("/mnt/mapr/user/mmoldovan/phosphosites/data/human_oggs_no_parals_mammalia_phylip/", "/mnt/mapr/user/mmoldovan/phosphosites/data/mammalia_species_for_paml.nwk", "/mnt/mapr/user/mmoldovan/phosphosites/data/human_oggs_no_parals_mammalia_trees/")



