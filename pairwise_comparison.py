"""
pairwise_comparison

analyzes BLAST results, builds pairwise alignments and prints out a 
list of conservative phosphorilated residues.

Input:
1. BLAST results in the output format 7
2. Two FASTA files with proteins
3. Two TSV files with phosphorilation info
"""

from optparse import OptionParser
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MuscleCommandline
from os import system, remove


def blastparse(blastout_filename):
	homol_dict = dict()
	with open(blastout_filename) as blastout:
		for s in blastout:
			if s[0] != '#' and not id_1:
				s = s.strip().split()
				id_1 = s[0].split('|')[2]
				id_2 = s[1].split('|')[2]
				homol_dict[id_1] = id_2
			if s[0] == '#':
				id_1 = ""
				id_2 = ""
	return homol_dict


def build_seqdict(fafile_name):
	seqdict = dict()
	with open(fafile_name) as fafile:
		for record in SeqIO.parse(fafile, "fasta"):
			seqdict[record.id.split('|')[-1]] = record
	return seqdict


def build_psite_dict(psite_filename):
	psite_dict = dict()
	with open(psite_filename) as psites:
		for s in psites:
			s = s.strip().split()
			if s[0] == "UniProtID":
				continue
			if psite_dict.get(s[0]):
				psite_dict[s[0]][eval(s[2])] = s
			else:
				psite_dict[s[0]] = dict()
				psite_dict[s[0]][eval(s[2])] = s
	return psite_dict


def seq_align(record1, record2):
	init_fasta = open("file_for_al_temp.fa", 'w')
	SeqIO.write(record1, init_fasta, "fasta")
	SeqIO.write(record2, init_fasta, "fasta")
	init_fasta.close()
#	cline = MuscleCommandline(input="file_for_al_temp.fa", out="file_al_temp.fa")
	cline = "muscle -in {} -out {} -quiet".format("file_for_al_temp.fa", "file_al_temp.fa")
	system(cline)


def align_work(psite_dict_1, psite_dict_2, uniprot_id_1, uniprot_id_2, n):
	alignment = AlignIO.read("file_al_temp.fa", "fasta")
	seq1 = ""
	seq2 = ""
	for record in alignment:
#		print record.id
#		print record.seq
#		print psite_dict_1[uniprot_id_1]
#		print psite_dict_2[uniprot_id_2]
#		print record.id
		if record.id.split('|')[-1] == uniprot_id_1:
			for c in record.seq:
				seq1 += c
		elif record.id.split('|')[-1] == uniprot_id_2:
			for c in record.seq:
				seq2 += c
	l = len(seq1)
	index_1 = 0
	index_2 = 0
#	print l
#	print uniprot_id_1
#	print seq1
#	print uniprot_id_2
#	print seq2
#	print psite_dict_1[uniprot_id_1]
#	print psite_dict_2[uniprot_id_2]
	for i in range(l):
		printed = []
		for j in range(index_1 - n, index_1 + n + 1):
#			if psite_dict_2[uniprot_id_2].get(index_2):
#				print "hello", index_2, '1'
			if psite_dict_1[uniprot_id_1].get(j) and psite_dict_2[uniprot_id_2].get(index_2):
				outstr = '\t'.join(psite_dict_1[uniprot_id_1][j]) + '\t' + '\t'.join(psite_dict_2[uniprot_id_2][index_2])
				if outstr not in printed:
					print outstr
					printed.append(outstr)


		if seq1[i] != '-':
			index_1 += 1
		if seq2[i] != '-':
			index_2 += 1
#		print index_1, index_2, seq1[i], seq2[i]


def main(blastout_filename, fafile_name_1, fafile_name_2, psite_filename_1, psite_filename_2, radius):
	seqdict_1 = build_seqdict(fafile_name_1)
	seqdict_2 = build_seqdict(fafile_name_2)
	psite_dict_1 = build_psite_dict(psite_filename_1)
	psite_dict_2 = build_psite_dict(psite_filename_2)
	homol_dict = blastparse(blastout_filename)
	for uniprot_id in homol_dict.keys():
		try:
			seq_align(seqdict_1[uniprot_id], seqdict_2[homol_dict[uniprot_id]])
			align_work(psite_dict_1, psite_dict_2, uniprot_id, homol_dict[uniprot_id], radius)
			remove("file_for_al_temp.fa")
			remove("file_al_temp.fa")
		except KeyError:
			pass


#		break

parser = OptionParser()
parser.add_option("-l", "--blastout_filename", help="BLAST results in the output format 7")
parser.add_option("-a", "--fafile_name_1", help="FASTA files with proteins 1 (BLAST query)")
parser.add_option("-b", "--fafile_name_2", help="FASTA files with proteins 2 (BLAST target)")
parser.add_option("-1", "--psite_filename_1", help="TSV files with phosphorilation info 1")
parser.add_option("-2", "--psite_filename_2", help="TSV files with phosphorilation info 2")
parser.add_option("-n", "--radius", help="Radius on the sequence for homologous p-sites")
opt, args = parser.parse_args()

main(opt.blastout_filename, opt.fafile_name_1, opt.fafile_name_2, opt.psite_filename_1, opt.psite_filename_2, eval(opt.radius))



