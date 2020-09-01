"""
prepare_als_and_trees

Prepare alignments and trees for the search for positive selection

Input:
1. Protein sequence directory
2. Nucleotide sequence alignment directory
3. Table with species IDs
4. Tree directory. The leaves must be named by species IDs

Output:
1. Directory with renamed nucleotide sequences
2. Directory with trimmed trees
"""

from os import listdir, mkdir
from ete3 import Tree
from optparse import OptionParser


def protID_spec_dict_build(protseq_dir):
	prot_id_spec_dict = dict()
	for filename in listdir(protseq_dir):
		with open(protseq_dir + filename) as inhandle:
			for s in inhandle:
				if s.startswith('>'):
					prot_id = s.split()[0][1:]
					spec = s.split('[')[-1].split(']')[0]
					prot_id_spec_dict[prot_id] = '_'.join(spec.split())
	return prot_id_spec_dict


def specID_dict_build(spec_id_file):
	spec_id_dict = dict()
	with open(spec_id_file) as inhandle:
		for s in inhandle:
			s = s.strip().split()
			spec_id_dict[s[0]] = s[1]
	return spec_id_dict


def protID_specID_dict_build(prot_id_spec_dict, spec_id_dict):
	prot_id_spec_id_dict = dict()
	for k in prot_id_spec_dict.keys():
		prot_id_spec_id_dict[k] = spec_id_dict[prot_id_spec_dict[k]]
	return prot_id_spec_id_dict


def rename_nucl_seqs(prot_id_spec_id_dict, nucl_file, outfile_name):
	spec_list = []
	with open(nucl_file) as inhandle, open(outfile_name, 'w') as outhandle:
		for s in inhandle:
			if s.startswith('>'):
				prot_id = s.strip()[1:]
				spec_id = prot_id_spec_id_dict[prot_id]
				outhandle.write('>' + spec_id + '\n')
				spec_list.append(spec_id)
			else:
				outhandle.write(s)
	return spec_list

def tree_prune(spec_list, tree_file, outfile_name):
	with open(tree_file) as t:
		for s in t:
			if s.startswith('('):
				tree_str = s
	with open("tmp.nwk", 'w') as tmpfile:
		tmpfile.write(s)
	t = Tree("tmp.nwk")
#	q = []
#	for treenode in spec_list:
#		print(t.get_leaves_by_name(treenode))
#		q.append(t.get_leaves_by_name(treenode))
#	a = []
#	for node in t.traverse("postorder"):
#		if node.is_leaf():
#			a.append(node.name)
#	print(a)
	t.prune(spec_list)
	t.write(format=1, outfile=outfile_name)


def main(protseq_dir, spec_id_file, nucl_al_dir, treefile_dir, newdir_name_nucl, newdir_name_tree):
	try:
		mkdir(newdir_name_nucl)
	except:
		pass
	try:
		mkdir(newdir_name_tree)
	except:
		pass
	prot_id_spec_dict = protID_spec_dict_build(protseq_dir)
	spec_id_dict = specID_dict_build(spec_id_file)
	prot_id_spec_id_dict = protID_specID_dict_build(prot_id_spec_dict, spec_id_dict)
	for filename in listdir(nucl_al_dir):
		spec_list = rename_nucl_seqs(prot_id_spec_id_dict, nucl_al_dir + filename, newdir_name_nucl + filename)
		tree_name = filename.split('.')[0] + ".nwk"
		tree_prune(spec_list, treefile_dir + tree_name, newdir_name_tree + tree_name)


parser = OptionParser()
parser.add_option("-n", "--nucl_al_dir", help="Nucleotide sequence alignment directory")
parser.add_option("-p", "--protseq_dir", help="Protein sequence directory")
parser.add_option("-s", "--spec_id_file", help="Table with species IDs")
parser.add_option("-t", "--treefile_dir", help="Tree directory. The leaves must be named by species IDs")
parser.add_option("-1", "--newdir_name_nucl", help="Future directory with renamed nucleotide sequences")
parser.add_option("-2", "--newdir_name_tree", help="Future directory with trimmed trees")
opt, args = parser.parse_args()

main(opt.protseq_dir, 
	 opt.spec_id_file, 
	 opt.nucl_al_dir, 
	 opt.treefile_dir, 
	 opt.newdir_name_nucl, 
	 opt.newdir_name_tree)

