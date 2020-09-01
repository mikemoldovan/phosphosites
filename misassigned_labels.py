"""
misassigned_labels

analyzes BLAST results, builds pairwise alignments and prints out a 
list of conservative phosphorilated residues.

Input:
1. BLAST results in the output format 7
2. Two FASTA files with proteins
3. Two TSV files with phosphorilation info
"""

from optparse import OptionParser
from Bio import SeqIO, AlignIO
from os import listdir

def blastparse(blastout_filename):
	homol_dict_1 = dict()
	homol_dict_2 = dict()
	with open(blastout_filename) as blastout:
		for s in blastout:
			if s[0] != '#' and not id_1:
				s = s.strip().split()
				id_1 = s[0].split('|')[2]
				id_2 = s[1].split('|')[2]
				homol_dict_1[id_1] = id_2
				homol_dict_2[id_2] = id_1
			if s[0] == '#':
				id_1 = ""
				id_2 = ""
	return homol_dict_1, homol_dict_2


def build_psite_dict(psite_filename):
	psite_dict = dict()
	with open(psite_filename) as psites:
		for s in psites:
			s = s.strip().split()
			if s[0] == "UniProtID":
				continue
			if not psite_dict.get(s[0]):
				psite_dict[s[0]] = dict()
			psite_dict[s[0]][eval(s[2])] = s
	return psite_dict


def align_work(align_file, psite_dict_1, psite_dict_2, uniprot_id_1, uniprot_id_2, n):
	alignment = AlignIO.read(align_file, "fasta")
	seq1 = ""
	seq2 = ""

	for record in alignment:
		if record.id.split('|')[-1] == uniprot_id_1:
			for c in record.seq:
				seq1 += c
		elif record.id.split('|')[-1] == uniprot_id_2:
			for c in record.seq:
				seq2 += c

	l = len(seq1)
	index_1 = 0
	index_2 = 0

	predictions = {"1->2_real" : 0,
				   "1->2_false" : 0,
				   "2->1_real" : 0,
				   "2->1_false" : 0}

	processed_inds_1 = [-1]
	processed_inds_2 = [-1]

	for i in range(l):
		for j in range(index_1 - n, index_1 + n + 1):
			if seq1[i] == '-' or seq2[i] == '-':
				continue

			if index_1 == processed_inds_1[-1]:
				pass
			elif psite_dict_1.get(uniprot_id_1):
				if psite_dict_1[uniprot_id_1].get(j) and \
				   psite_dict_2.get(uniprot_id_2) and \
				   psite_dict_2[uniprot_id_2].get(index_2):

					predictions["1->2_real"] += 1
				elif psite_dict_1[uniprot_id_1].get(j):
					predictions["1->2_false"] += 1

			if index_2 == processed_inds_2[-1]:
				pass
			elif psite_dict_2.get(uniprot_id_2):
				if psite_dict_1.get(uniprot_id_1) and \
				   psite_dict_1[uniprot_id_1].get(j) and \
				   psite_dict_2[uniprot_id_2].get(index_2):

					predictions["2->1_real"] += 1
				elif psite_dict_2[uniprot_id_2].get(j):
					predictions["2->1_false"] += 1

		processed_inds_1.append(index_1)
		processed_inds_2.append(index_2)

		if seq1[i] != '-':
			index_1 += 1
		if seq2[i] != '-':
			index_2 += 1

	return predictions


def sum_dicts(dict1, dict2):
	for k in dict1.keys():
		dict1[k] += dict2[k]
	return dict1


def main(blastout_filename, aldir_name, psite_filename_1, psite_filename_2, radius):
	predictions = {"1->2_real" : 0.,
				   "1->2_false" : 0.,
				   "2->1_real" : 0.,
				   "2->1_false" : 0.}

	psite_dict_1 = build_psite_dict(psite_filename_1)
	psite_dict_2 = build_psite_dict(psite_filename_2)
	homol_dict_1, homol_dict_2 = blastparse(blastout_filename)

	for align_file in listdir(aldir_name):
		id_1 = '.'.join(align_file.split('.')[:-1])
#		if not homol_dict_1.get(id_1):
#			continue
		id_2 = homol_dict_1[id_1]
#		print(id_1, id_2)
		align_file_path = aldir_name + align_file
		p = align_work(align_file_path, psite_dict_1, psite_dict_2, id_1, id_2, radius)
		predictions = sum_dicts(predictions, p)

	p1 = predictions["1->2_real"]/(predictions["1->2_real"] + predictions["1->2_false"])
	p2 = predictions["2->1_real"]/(predictions["2->1_real"] + predictions["2->1_false"])
	print("1->2\t{}".format(p1))
	print("2->1\t{}".format(p2))

	print(predictions)


parser = OptionParser()
parser.add_option("-l", "--blastout_filename", help="BLAST results in the output format 7")
parser.add_option("-a", "--aldir_name", help="directory with alignment data")
parser.add_option("-1", "--psite_filename_1", help="TSV files with phosphorilation info 1")
parser.add_option("-2", "--psite_filename_2", help="TSV files with phosphorilation info 2")
parser.add_option("-n", "--radius", help="Radius on the sequence for homologous p-sites")
opt, args = parser.parse_args()

main(opt.blastout_filename, opt.aldir_name, opt.psite_filename_1, opt.psite_filename_2, eval(opt.radius))

