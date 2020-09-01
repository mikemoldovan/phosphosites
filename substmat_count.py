"""
substmat_count

Count substitution matrices distinguishing phosphorilated
and not phosphorilated sites

Input:
1. Directory with paml_run outputs
2. File containing a list of phosphosites
3. Identifier of species with phosphosite data

Output:
substitution matrix
"""


from optparse import OptionParser
from os import listdir
from Bio import SeqIO
from ete3 import Tree


def make_psite_dict(psite_filename):
	psite_dict = dict()
	with open(psite_filename) as psite_file:
		for s in psite_file:
			s = s.strip().split()
			if s[0] == "UniProtID":
				continue
			if psite_dict.get(s[0]):
				psite_dict[s[0]][s[2]] = s[3]
			else:
				psite_dict[s[0]] = dict()
				psite_dict[s[0]][s[2]] = s[3]
	return psite_dict


def make_seq_dict(fafile_name):
	seq_dict = dict()
	with open(fafile_name) as fafile:
		for record in SeqIO.parse(fafile, "fasta"):
			seq = ""
			for c in record.seq:
				seq += c
			seq_dict[record.id] = seq
	return seq_dict


#retunrns psite coordinate dict (coordinate start at 0)
def infer_psites_in_al(seq_dict, psite_dict, spec_id, gene_id):
	psites_dict = dict()
	for k in seq_dict.keys():
		if k.split('_')[-1] == spec_id:
			psite_coord_dict = dict()
			count = 0
			count2 = 0
			for c in seq_dict[k]:
				if c != '-':
					count += 1
				if psite_dict[gene_id].get(str(count)):
					psites_dict[count2] = True
				count2 += 1
#			print list(sorted(psites_dict.keys()))
#			print map(eval, list(sorted(psite_dict[gene_id].keys())))
#			print seq_dict[k]
			return psites_dict
	return dict()


#seq1 -- ancestral, seq2 -- descendal
def compare_align_seqs(seq1, seq2, seqlen, psites_dict, subst_mat):
	for i in range(seqlen):
		c1 = seq1[i]
		c2 = seq2[i]
		if psites_dict.get(i):
			if c1 in ('T', 'Y', 'S'):
				c1 = 'p' + c1
			if c2 in ('T', 'Y', 'S'):
				c2 = 'p' + c2
		if c1 == '-' or c1 == 'X':
			continue
		if c2 == '-' or c2 == 'X':
			continue
		if not subst_mat.get(c1):
			subst_mat[c1] = dict()
		if not subst_mat[c1].get(c2):
			subst_mat[c1][c2] = 0
		subst_mat[c1][c2] += 1


def align_analyze(gene_id, paml_run_outdir, spec_id, psite_dict, subst_mat):
	with open(paml_run_outdir + gene_id + ".nwk") as treefile:
		for s in treefile:
			s = s.strip()
			if s.startswith('('):
				tree_str = s
	seq_dict = make_seq_dict(paml_run_outdir + gene_id + ".fa")
	try:
		psites_dict = infer_psites_in_al(seq_dict, psite_dict, spec_id, gene_id)
	except KeyError:
		psites_dict = dict()
	t = Tree(tree_str, format=8)
	for node in t.traverse("postorder"):
		if node.is_leaf():
			continue
		seq1 = seq_dict[node.name]
		seqlen = len(seq1)
		for child in node.children:
			seq2 = seq_dict[child.name]
			compare_align_seqs(seq1, seq2, seqlen, psites_dict, subst_mat)


#column -- ancestor, row -- descendant
def print_subst_mat(subst_mat):
	colnames = sorted(list(subst_mat.keys()))
	colnames_l = len(colnames)
	rownames = []
	for k1 in subst_mat.keys():
		for k2 in subst_mat[k1].keys():
			rownames.append(k2)
	rownames = sorted(list(set(rownames)))
	rownames_l = len(rownames)
	print '\t'.join(['#'] + colnames)
	for c1 in colnames:
		row = [c1]
		for c2 in rownames:
			if subst_mat[c1].get(c2):
				row.append(subst_mat[c1][c2])
			else:
				row.append(0)
		print '\t'.join(map(str, row))


def main(paml_run_outdir, psite_filename, spec_id):
	psite_dict = make_psite_dict(psite_filename)
	subst_mat = dict()
	count = 0
	for filename in listdir(paml_run_outdir):
		if filename.endswith(".fa"):
			gene_id = filename[:-3]
			align_analyze(gene_id, paml_run_outdir, spec_id, psite_dict, subst_mat)
			count += 1
			if count % 1000 == 0:
				print subst_mat
				print_subst_mat(subst_mat)
	print subst_mat
	print_subst_mat(subst_mat)


parser = OptionParser()
parser.add_option("-p", "--paml_run_outdir", help="Directory with paml_run outputs")
parser.add_option("-f", "--psite_filename", help="File containing a list of phosphosites")
parser.add_option("-i", "--spec_id", help="Identifier of species with phosphosite data")
opt, args = parser.parse_args()


main(opt.paml_run_outdir, opt.psite_filename, opt.spec_id)

#main("/home/mmoldovan/phosphosites/data/human_oggs_paml_out/", "/home/mmoldovan/phosphosites/data/human.tsv", "4054")
#main("/home/mmoldovan/phosphosites/data/human_oggs_paml_out/", "/home/mmoldovan/phosphosites/data/human_minscore2.tsv", "4054")
#main("/home/mmoldovan/phosphosites/data/human_oggs_paml_out/", "/home/mmoldovan/phosphosites/pair_compar/mouse_human_0_2.tsv", "4054")
#main("/home/mmoldovan/phosphosites/data/human_oggs_paml_out/", "/home/mmoldovan/phosphosites/data/human_in_disordered.tsv", "4054")
#main("/home/mmoldovan/phosphosites/data/human_oggs_paml_out/", "/home/mmoldovan/phosphosites/data/human_in_ordered.tsv", "4054")
#"/home/mmoldovan/phosphosites/data/human_minscore2.tsv"
#"/home/mmoldovan/phosphosites/pair_compar/mouse_human_0.tsv", "2597"

human_in_disordered.tsv
human_in_ordered.tsv
mouse_human_0_2.tsv

MTGYTPDEKLRLDQLRALRRQWLKDQELSPREPVLPPEKAWPFSGFWHRFLQKDSAWRRFAYKAYGSGMFVFVNFLIPAWIVHYYVKYHVETRPYGIVETKRKIFPGDTILETGEVIPPMKEHDTHNH
