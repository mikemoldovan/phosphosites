"""
phosphosite_retention

analyzes BLAST results, builds pairwise alignments and, 
for a given aminoacid pair (say, S and T), prints out 
the following values: 
#[pS - pS]
#[pS - S]
#[pS - pT]
#[pS - T]

and the p-value of contingency matrix:

#[pS - pS]				#[pS - pT]
#[pS - pS] + #[pS - S]	#[pS - pT] + #[pS - T]

Input:
1. BLAST results in the output format 7
2. Two FASTA files with proteins
3. Two TSV files with phosphorilation info
4. Two aminoacids
"""

from optparse import OptionParser
from Bio import SeqIO, AlignIO
from os import system, remove, mkdir, path
#from pathlib import Path
#from scipy import stats


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
				psite_dict[s[0]][eval(s[2]) - 1] = s
			else:
				psite_dict[s[0]] = dict()
				psite_dict[s[0]][eval(s[2]) - 1] = s
	return psite_dict


def seq_align(record1, record2, uniprot_id, alignment_dir):
	if path.isfile(alignment_dir + uniprot_id + ".fa"):
		return alignment_dir + uniprot_id + ".fa"
	init_fasta = open("file_for_al_temp.fa", 'w')
	SeqIO.write(record1, init_fasta, "fasta")
	SeqIO.write(record2, init_fasta, "fasta")
	init_fasta.close()
	cline = "muscle -in {} -out {} -quiet".format("file_for_al_temp.fa", alignment_dir + uniprot_id + ".fa")
	system(cline)
	return alignment_dir + uniprot_id + ".fa"


def align_work(psite_dict_1, psite_dict_2, uniprot_id_1, uniprot_id_2, aas, outdict, alname):
	alignment = AlignIO.read(alname, "fasta")
	seq1 = ""
	seq2 = ""
	for record in alignment:
		if record.id.split('|')[-1] == uniprot_id_1:
			for c in record.seq:
				seq1 += c
		elif record.id.split('|')[-1] == uniprot_id_2:
			for c in record.seq:
				seq2 += c
	if len(seq1) != len(seq2):
		print uniprot_id_1, uniprot_id_2
		print "seq1", seq1
		print "seq2", seq2
	l = len(seq1)
	index_1 = 0
	index_2 = 0
	ind_1_list = dict()
	ind_2_list = dict()
	for i in range(l):
		c_get1 = psite_dict_1.get(uniprot_id_1)
		c_get2 = psite_dict_2.get(uniprot_id_2)
		if c_get1:
			ispsite_1 = psite_dict_1[uniprot_id_1].get(index_1)
		else:
			ispsite_1 = False
		if c_get2:
			ispsite_2 = psite_dict_2[uniprot_id_2].get(index_2)
		else:
			ispsite_2 = False
		if c_get1 and c_get2 and ispsite_1 and ispsite_2:
			if ind_1_list.get(index_1):
				continue
			if ind_2_list.get(index_2):
				continue
			ind_1_list[index_1] = True
			ind_2_list[index_2] = True
			if seq1[i] == aas[0] and seq2[i] == aas[0]:
				outdict["#pX1|pX1"] += 1
			elif seq1[i] == aas[0] and seq2[i] == aas[1]:
				outdict["#pX1|pX2"] += 1
			elif seq1[i] == aas[1] and seq2[i] == aas[0]:
				outdict["#pX1|pX2"] += 1
			elif seq1[i] == aas[1] and seq2[i] == aas[1]:
				outdict["#pX2|pX2"] += 1
		
		elif c_get1 and ispsite_1:
			if ind_1_list.get(index_1):
				continue
			ind_1_list[index_1] = True
			if seq1[i] == aas[0] and seq2[i] == aas[0]:
				outdict["#pX1|X1"] += 1
			elif seq1[i] == aas[0] and seq2[i] == aas[1]:
				outdict["#pX1|X2"] += 1
			elif seq1[i] == aas[1] and seq2[i] == aas[0]:
				outdict["#X1|pX2"] += 1
			elif seq1[i] == aas[1] and seq2[i] == aas[1]:
				outdict["#pX2|X2"] += 1

		elif c_get2 and ispsite_2:
			if ind_2_list.get(index_2):
				continue
			ind_2_list[index_2] = True
			if seq1[i] == aas[0] and seq2[i] == aas[0]:
				outdict["#pX1|X1"] += 1
			elif seq1[i] == aas[0] and seq2[i] == aas[1]:
				outdict["#X1|pX2"] += 1
			elif seq1[i] == aas[1] and seq2[i] == aas[0]:
				outdict["#pX1|X2"] += 1
			elif seq1[i] == aas[1] and seq2[i] == aas[1]:
				outdict["#pX2|X2"] += 1
		if len(seq1) != len(seq2):
			print uniprot_id_1, uniprot_id_2
			print "seq1", seq1
			print "seq2", seq2

		if seq1[i] != '-':
			index_1 += 1
		if seq2[i] != '-':
			index_2 += 1		


def main(blastout_filename, fafile_name_1, fafile_name_2, psite_filename_1, psite_filename_2, aas, alignment_dir):
	try:
		mkdir(alignment_dir)
	except:
		pass
	seqdict_1 = build_seqdict(fafile_name_1)
	seqdict_2 = build_seqdict(fafile_name_2)
	psite_dict_1 = build_psite_dict(psite_filename_1)
	psite_dict_2 = build_psite_dict(psite_filename_2)
	homol_dict = blastparse(blastout_filename)
	outdict = {"#pX1|pX1" : 0, #
			   "#pX1|X1" : 0,
			   "#pX1|pX2" : 0, #
			   "#pX1|X2" : 0,
			   "#X1|pX2" : 0,
			   "#pX2|pX2" : 0,
			   "#pX2|X2" : 0} #
	count = 0
	l = len(list(homol_dict.keys()))
	for uniprot_id in homol_dict.keys():
#		if count > 100:
#			break
		count += 1
#		try:
		if not seqdict_1.get(uniprot_id):
			continue
#			print '1', uniprot_id
		if not seqdict_2.get(homol_dict[uniprot_id]):
			continue
#			print '2', homol_dict[uniprot_id]
		if seqdict_1.get(uniprot_id) and seqdict_2.get(homol_dict[uniprot_id]):
			alname = seq_align(seqdict_1[uniprot_id], seqdict_2[homol_dict[uniprot_id]], uniprot_id, alignment_dir)
#			print l, count, uniprot_id, homol_dict[uniprot_id], seqdict_1[uniprot_id].id, seqdict_2[homol_dict[uniprot_id]].id
			try:
				align_work(psite_dict_1, psite_dict_2, uniprot_id, homol_dict[uniprot_id], aas, outdict, alname)
			except ValueError:
				print outdict
			try:
				remove("file_for_al_temp.fa")
				remove("file_al_temp.fa")
			except:
				pass
#		except KeyError:
#			pass

#	statfunc = stats.chi2_contingency
#	cont_table1 = [[outdict["#pX1|pX1"], outdict["#pX1|pX1"] + outdict["#pX1|pX2"]], [outdict["#pX1|pX2"], outdict["#pX1|pX2"] + outdict["#pX1|X2"]]]
#	r1 = float(cont_table1[0][0])/cont_table1[0][1]
#	r2 = float(cont_table1[1][0])/cont_table1[1][1]
#	g, pvalue, dof, expctd = statfunc(cont_table1)
#	print "{}\t{}\t{}\t{}".format(aas[0], r1, r2, pvalue)
	print outdict

#	cont_table2 = [[outdict["#pX2|pX2"], outdict["#pX2|pX2"] + outdict["#pX1|pX2"]], [outdict["#pX1|pX2"], outdict["#pX1|pX2"] + outdict["#pX1|X2"]]]
#	r1 = float(cont_table1[0][0])/cont_table1[0][1]
#	r2 = float(cont_table1[1][0])/cont_table1[1][1]
#	g, pvalue, dof, expctd = statfunc(cont_table1)


#		break

parser = OptionParser()
parser.add_option("-l", "--blastout_filename", help="BLAST results in the output format 7")
parser.add_option("-a", "--fafile_name_1", help="FASTA files with proteins 1 (BLAST query)")
parser.add_option("-b", "--fafile_name_2", help="FASTA files with proteins 2 (BLAST target)")
parser.add_option("-1", "--psite_filename_1", help="TSV files with phosphorilation info 1")
parser.add_option("-2", "--psite_filename_2", help="TSV files with phosphorilation info 2")
parser.add_option("-n", "--aas", help="Aminoacid pair")
parser.add_option("-f", "--alignment_dir", help="Directory with pairwise alignments")
opt, args = parser.parse_args()

main(opt.blastout_filename, opt.fafile_name_1, opt.fafile_name_2, opt.psite_filename_1, opt.psite_filename_2, opt.aas, opt.alignment_dir)



"""
python ../scripts/phosphosite_retention.py -l mouse_big.txt -a ../data/proteomes/mouse_big_proteome.fa -b ../data/proteomes/human_big_proteome.fa -1 ../data/sitelists_filtered/mouse_iptmnet.tsv -2 ../data/sitelists_filtered/human_iptmnet.tsv -n ST -f human_mouse_als/ > mouse_human_ST.txt
python ../scripts/phosphosite_retention.py -l mouse_big.txt -a ../data/proteomes/mouse_big_proteome.fa -b ../data/proteomes/human_big_proteome.fa -1 ../data/sitelists_filtered/mouse_iptmnet.tsv -2 ../data/sitelists_filtered/human_iptmnet.tsv -n SY -f human_mouse_als/ > mouse_human_SY.txt
python ../scripts/phosphosite_retention.py -l mouse_big.txt -a ../data/proteomes/mouse_big_proteome.fa -b ../data/proteomes/human_big_proteome.fa -1 ../data/sitelists_filtered/mouse_iptmnet.tsv -2 ../data/sitelists_filtered/human_iptmnet.tsv -n TY -f human_mouse_als/ > mouse_human_TY.txt
/home/mmoldovan/octopus_rna_editing/data/orfs_bim_plus.fa
/home/mmoldovan/octopus_rna_editing/data/orfs_oct_plus.fa
/home/mmoldovan/octopus_rna_editing/data/orfs_squ_plus.fa
/home/mmoldovan/octopus_rna_editing/data/orfs_sep_plus.fa
/home/mmoldovan/octopus_rna_editing/data/orfs_nau_plus.fa
/home/mmoldovan/octopus_rna_editing/data/orfs_apl_plus.fa

/home/mmoldovan/octopus_rna_editing/data/editing_sites_bim_plus.txt
/home/mmoldovan/octopus_rna_editing/data/editing_sites_oct_plus.txt
/home/mmoldovan/octopus_rna_editing/data/editing_sites_squ_plus.txt
/home/mmoldovan/octopus_rna_editing/data/editing_sites_sep_plus.txt
/home/mmoldovan/octopus_rna_editing/data/editing_sites_nau_plus.txt
/home/mmoldovan/octopus_rna_editing/data/editing_sites_apl_plus.txt
"""
