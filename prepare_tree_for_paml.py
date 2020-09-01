"""
prepare_tree_for_paml

Replace leaf names in a tree with numerical identifiers

Input:
Tree in Newick format

Output:
1. Tree with leaf names substituted with identifiers
2. Table with identifiers and species names
"""


from optparse import OptionParser
from Bio import Phylo


def main(nw_tree_filename):
	tree = Phylo.read(nw_tree_filename, "newick")
	count = 0
	for leaf in tree.get_terminals():
		print leaf.name + '\t' + str(count)
		leaf.name = str(count)
		count += 1
	new_tree_name = '.'.join(nw_tree_filename.split('.')[:-1]) + "_for_paml.nwk"
	Phylo.write(tree, new_tree_name, "newick")


parser = OptionParser()
parser.add_option("-i", "--nw_tree_filename", help="Tree in Newick format")
opt, args = parser.parse_args()


main(opt.nw_tree_filename)