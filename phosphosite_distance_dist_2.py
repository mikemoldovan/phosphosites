"""
phosphosite_distance_dist_2

Count distributions of normalized phosphosite densities and distances
between phosphosites for different disordered regions

Input:
1. Table with phosphosites
2. Table with disordered regions
3. File with proteome
4. Number of permutations 
5. Aminoacids for which phosphosites to be counted (default 'STY')
6. Upper threshold for phosphosite coverage (default 0.5)

Output:
1. Table with normalized distribution of phosphosite densities
2. Table with permutations for [1]
3. Table with distribution of distances between phosphosites
4. Table with permutations for [2]
"""


import random
from optparse import OptionParser
from Bio import SeqIO


def make_psite_dict(psite_filename, seq_dict, aacid_list):
	psite_dict = dict()
	with open(psite_filename) as psite_file:
		for s in psite_file:
			s = s.strip().split()
			if s[0] == "UniProtID":
				continue
			if seq_dict.get(s[0]):
				seq = seq_dict[s[0]]
				if len(seq) <= eval(s[2]) - 1:
					continue
				if seq[eval(s[2]) - 1] not in aacid_list:
					continue
			else:
				continue
			if psite_dict.get(s[0]):
				psite_dict[s[0]][eval(s[2]) - 1] = s[3]
			else:
				psite_dict[s[0]] = dict()
				psite_dict[s[0]][eval(s[2]) - 1] = s[3]
	return psite_dict


def make_disord_reg_dict(disord_reg_filename):
	disord_reg_dict = dict()
	with open(disord_reg_filename) as drf:
		for s in drf:
			s = s.strip().split()
			if not s:
				continue
			disord_regs = []
			if len(s) < 2:
				continue
			for crd_pair in s[1].split('|'):
				crd_pair = crd_pair.split('-')
				disord_regs.append([eval(crd_pair[0]) - 1, eval(crd_pair[1])])
			disord_reg_dict[s[0]] = disord_regs
	return disord_reg_dict


def make_seq_dict(proteome_fasta):
	seq_dict = dict()
	with open(proteome_fasta) as infasta:
		for record in SeqIO.parse(infasta, "fasta"):
			seq = ""
			for c in record.seq:
				seq += c
			seq_dict[record.id.split('|')[2]] = seq
	return seq_dict


def make_pos_lists(seq, psite_dict_seq, disord_reg_list, aacid_list):
	pos_list = []
	psite_pos_list = []
	curr_disord_reg = 0
	dr_pos_list = []
	dr_psite_pos_list = []
	len_disord_list = len(disord_reg_list)
	for i in range(len(seq)):
		if i >= disord_reg_list[curr_disord_reg][1]:
			curr_disord_reg += 1
			if curr_disord_reg == len_disord_list:
				break
			if dr_pos_list:
				pos_list.append(dr_pos_list)
				psite_pos_list.append(dr_psite_pos_list)
				dr_pos_list = []
				dr_psite_pos_list = []

		if i >= disord_reg_list[curr_disord_reg][0] and i < disord_reg_list[curr_disord_reg][1]:
			if seq[i] in aacid_list:
				dr_pos_list.append(i)
				if psite_dict_seq.get(i):
					dr_psite_pos_list.append(i)
	if dr_pos_list:
		pos_list.append(dr_pos_list)
		psite_pos_list.append(dr_psite_pos_list)

	return pos_list, psite_pos_list


def shuffle_poss(pos_list, psite_pos_list, permut_num):
	shuff_pos_lists = []
	shuff_psite_pos_lists = []
	l = len(pos_list)
	positions = []
#	psite_positions = []
	psite_num = 0
	for i in range(l):
		for pos in pos_list[i]:
			positions.append(pos)
			if pos in psite_pos_list[i]:
				psite_num += 1
	for i in range(permut_num):
		rand_psites = sorted(random.sample(positions, psite_num))
		shuff_pos_list = []
		for i in range(l):
			dr_pos_list = []
			for j in range(psite_num):
				if rand_psites[j] in pos_list[i]:
					dr_pos_list.append(rand_psites[j])
			shuff_pos_list.append(dr_pos_list)
		shuff_pos_lists.append(shuff_pos_list)
	return shuff_pos_lists


def stat_count(pos_list, psite_pos_list, site_number):
	S1_list = []
	S2_list = []
	l = len(pos_list)
	for i in range(l):
		l_psite = float(len(psite_pos_list[i]))
#		S1 = (l_psite/len(pos_list[i]))/float(site_number)
		if len(pos_list[i]) != 0:
			S1 = l_psite/len(pos_list[i])
			S1_list.append(S1)
		if l_psite > 1:
			for j in range(1, int(l_psite)):
				try:
					S2 = pos_list[i].index(psite_pos_list[i][j]) - pos_list[i].index(psite_pos_list[i][j - 1])
					S2_list.append(S2)
				except:
					pass
	return S1_list, S2_list


def site_num_count(pos_list):
	site_num = 0
	for arr in pos_list:
		site_num += len(arr)
	return site_num


def main(psite_filename,
		 disord_reg_filename,
		 proteome_fasta,
		 permut_num,
		 aacid_list,
		 upper_thr,
		 job_id):
	aacid_list = list(aacid_list)
	disord_reg_dict = make_disord_reg_dict(disord_reg_filename)
	seq_dict = make_seq_dict(proteome_fasta)
	psite_dict = make_psite_dict(psite_filename, seq_dict, aacid_list)

	S1_outfile = open("{}S1_stats.txt".format(job_id), 'w')
	S2_outfile = open("{}S2_stats.txt".format(job_id), 'w')
	S1_permuts = open("{}S1_permuts.txt".format(job_id), 'w')
	S2_permuts = open("{}S2_permuts.txt".format(job_id), 'w')

	for uniprot_id in seq_dict.keys():
		try:
			i = seq_dict[uniprot_id]
			i = psite_dict[uniprot_id]
			i = disord_reg_dict[uniprot_id]
		except:
			continue

		pos_list, psite_pos_list = make_pos_lists(seq_dict[uniprot_id],
												  psite_dict[uniprot_id], 
												  disord_reg_dict[uniprot_id], 
												  aacid_list)
		psite_pos_site_num = site_num_count(psite_pos_list)
		pos_site_num = site_num_count(pos_list)
		if pos_site_num == 0:
			continue
		if float(psite_pos_site_num)/pos_site_num > upper_thr:
			continue
		if psite_pos_site_num < 2:
			continue
		S1_list, S2_list = stat_count(pos_list, psite_pos_list, pos_site_num)
		for i in S1_list:
			S1_outfile.write(str(i) + '\n')
		for i in S2_list:
			S2_outfile.write(str(i) + '\n')
		shuff_pos_lists = shuffle_poss(pos_list, psite_pos_list, permut_num)
		for shuff_pos_list in shuff_pos_lists:
			S1_list, S2_list = stat_count(pos_list, shuff_pos_list, pos_site_num)
#			print S1_list, S2_list
			for i in S1_list:
				S1_permuts.write(str(i) + '\n')
			for i in S2_list:
				S2_permuts.write(str(i) + '\n')

	S1_outfile.close()
	S2_outfile.close()
	S1_permuts.close()
	S2_permuts.close()

	print "{}S1_stats.txt".format(job_id) 
	print "{}S2_stats.txt".format(job_id)
	print "{}S1_permuts.txt".format(job_id)
	print "{}S2_permuts.txt".format(job_id) 


parser = OptionParser()
parser.add_option("-p", "--psite_filename", help="Table with phosphosites")
parser.add_option("-r", "--disord_reg_filename", help="Table with disordered regions")
parser.add_option("-f", "--proteome_fasta", help="File with proteome")
parser.add_option("-n", "--permut_num", help="Number of permutations (default='2')", default='2')
parser.add_option("-a", "--aacid_list", help="Aminoacids for which phosphosites to be counted (default 'STY')", default="STY")
parser.add_option("-u", "--upper_thr", help="Upper threshold for phosphosite coverage (default 0.5)", default='0.5')
parser.add_option("-j", "--job_id", help="Job identifier (default '')", default='')
opt, args = parser.parse_args()


main(opt.psite_filename,
	 opt.disord_reg_filename,
	 opt.proteome_fasta,
	 eval(opt.permut_num),
	 opt.aacid_list,
	 eval(opt.upper_thr),
	 opt.job_id)

