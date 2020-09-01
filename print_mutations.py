"""
print_mutations

Infer mutations from the PAML output and print the mutation list

Prints file in format:
#gene_id
pos_in_al	nodename_ancetral	nodename_descandal	pos_in_seq_ancestral	pos_in_seq_descendal	let_ancestral	let_descendal

Positions start with 0

input:
PAML output dir
"""

from optparse import OptionParser
from os import listdir
from Bio import SeqIO
from ete3 import Tree
from disord_reg_adlib import *


def make_seq_dict(fafile_name, selected_spec_id):
	seq_dict = dict()
	selected_seq_dict = dict()
	with open(fafile_name) as fafile:
		for record in SeqIO.parse(fafile, "fasta"):
			seq = ""
			for c in record.seq:
				seq += c
			seq_dict[record.id] = seq
			id_list = record.id.split('_')
			if len(id_list) == 2:
				if id_list[1] == selected_spec_id:
					count_seq = 0
#					print(seq)
					l = len(seq)
#					print(seq)
					al_index = 0
					for i in range(l):
						while al_index < l and seq[al_index] == '-':
							selected_seq_dict[al_index] = -1
							al_index += 1
						selected_seq_dict[al_index] = count_seq
						count_seq += 1
						al_index += 1
#	if selected_seq_dict:
#		print(selected_seq_dict)
	return seq_dict, selected_seq_dict


#seq1 -- ancestral, seq2 -- descendal
def compare_align_seqs(seq1, seq2, nodename_ancetral, nodename_descandal, seqlen, selected_seq_dict, dr_list):
	pos_in_prot_seq1 = 0
	pos_in_prot_seq2 = 0
	for i in range(seqlen):
		c1 = seq1[i]
		c2 = seq2[i]
		if c1 == '-' or c1 == 'X':
			if c1 == 'X':
				pos_in_prot_seq1 += 1
			continue
		if c2 == '-' or c2 == 'X':
			if c2 == 'X':
				pos_in_prot_seq2 += 1
			continue
		if c1 != c2:
#			print(selected_seq_dict)
			if selected_seq_dict.get(i):
				pos_in_selspec = selected_seq_dict[i]
			else:
				pos_in_selspec = 'N'
			if not selected_seq_dict.get(i):
				dr = 'N'
			elif selected_seq_dict[i] == -1:
				dr = '*'
			else:
				dr = pos_in_disord_reg(selected_seq_dict[i], dr_list)
				if not dr:
					dr = '*'
				else:
					dr = str(dr[0]) + '-' + str(dr[1])
			print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(str(i),
													nodename_ancetral,
													nodename_descandal,
													str(pos_in_prot_seq1),
													str(pos_in_prot_seq2),
													c1,
													c2,
													dr,
													pos_in_selspec))
		pos_in_prot_seq1 += 1
		pos_in_prot_seq2 += 1


def align_analyze(gene_id, paml_run_outdir, dr_list, selected_spec_id):
	print('#' + gene_id)
	with open(paml_run_outdir + gene_id + ".nwk") as treefile:
		for s in treefile:
			s = s.strip()
			if s.startswith('('):
				tree_str = s
	seq_dict, selected_seq_dict = make_seq_dict(paml_run_outdir + gene_id + ".fa", selected_spec_id)
#	selseq_al_drs = alseq_disord_crds(selected_seq, dr_list)
	t = Tree(tree_str, format=8)
	for node in t.iter_descendants():
		if node.is_leaf():
			continue
		seq1 = seq_dict[node.name]
		seqlen = len(seq1)
		for child in node.children:
			seq2 = seq_dict[child.name]
			compare_align_seqs(seq1, seq2, node.name, child.name, seqlen, selected_seq_dict, dr_list)


def main(paml_run_outdir, disord_reg_filename, selected_spec_id):
	dr_dict = build_dr_dict(disord_reg_filename)
	print("#pos_in_al\tnodename_ancetral\tnodename_descandal\tpos_in_seq_ancestral\tpos_in_seq_descendal\tlet_ancestral\tlet_descendal\tdisord_reg\tpos_in_ref_seq")
	for filename in listdir(paml_run_outdir):
		if filename.endswith(".fa"):
			gene_id = filename[:-3]
			align_analyze(gene_id, paml_run_outdir, dr_dict[gene_id], selected_spec_id)


parser = OptionParser()
parser.add_option("-i", "--paml_run_outdir", help="Directory with paml_run outputs")
parser.add_option("-d", "--disord_reg_filename", help="File with disordered regions")
parser.add_option("-s", "--selected_spec_id", help="ID of selected species (Human 4054, Mouse 2597)")
opt, args = parser.parse_args()

print("#disord_reg: N1-N2 for disordered region in reference species, * for ordered region,  N for no sequnce from reference species")
print("#tpos_in_ref_seq: N1 for position in the reference species sequence in alignment, -1 for gap, N for no sequnce from reference species")
main(opt.paml_run_outdir, opt.disord_reg_filename, opt.selected_spec_id)

#make_seq_dict(opt.paml_run_outdir + "ANXA6_HUMAN" + ".fa", opt.selected_spec_id)

#main("/home/mmoldovan/phosphosites/data/human_oggs_paml_out/", "/home/mmoldovan/phosphosites/data/human_oggs_mammalia_paml_mutations.txt")
#ANXA6_HUMAN	1-11|20-20|52-58|167-175|226-228|298-300|303-327|361-379|396-402|417-419|440-455|477-495|506-523|535-542|577-582|612-620|652-655|670-673

