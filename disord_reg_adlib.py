"""
disord_reg_adlib

Additional library helping to work with disordered regions
"""

from Bio import SeqIO

def build_dr_dict(dr_filename):
	dr_dict = dict()
	with open(dr_filename) as drfile:
		for dr_record in drfile:
			dr_record = dr_record.strip().split()
			dr_crds = []
			if len(dr_record) != 2:
				dr_dict[dr_record[0]] = dr_crds
				continue
			for crd_pair in dr_record[1].split('|'):
				crd_pair = crd_pair.split('-')
				dr_crds.append([eval(crd_pair[0]) - 1, eval(crd_pair[1])])
			dr_dict[dr_record[0]] = dr_crds
	return dr_dict


def make_psite_dict(psite_filename):
	psite_dict = dict()
	with open(psite_filename) as psite_file:
		for s in psite_file:
			s = s.strip().split()
			if s[0] == "UniProtID":
				continue
			if psite_dict.get(s[0]):
				psite_dict[s[0]][eval(s[2]) - 1] = s[3]
			else:
				psite_dict[s[0]] = dict()
				psite_dict[s[0]][eval(s[2]) - 1] = s[3]
	return psite_dict


def make_psite_dict2(psite_filename):
	psite_dict = dict()
	with open(psite_filename) as psite_file:
		for s in psite_file:
			s = s.strip().split()
			if s[0] == "UniProtID":
				continue
			if psite_dict.get(s[0]):
				psite_dict[s[0]][eval(s[2]) - 1] = s
			else:
				psite_dict[s[0]] = dict()
				psite_dict[s[0]][eval(s[2]) - 1] = s
	return psite_dict


#pos_num -- starts with 0
def pos_in_disord_reg(pos_num, dr_list):
	for crd_pair in dr_list:
		if pos_num < crd_pair[1] and pos_num >= crd_pair[0]:
			return crd_pair
	return None


def pos_in_al(alseq, seq_crd):
	seq_index = 0
	for al_index in range(len(alseq)):
		while alseq[al_index] == '-':
			al_index += 1
		if seq_crd == seq_index:
			return al_index
		seq_crd += 1


#alseq -- sequence with gaps
#dr_list -- list of disorderd regions for unaligned sequence
def alseq_disord_crds(alseq, dr_list):
	if dr_list == []:
		return []
	alseq_dr_crds = []
	seq_index = 0
	curr_dr_num = 0
	l = len(alseq)
	for al_index in range(l):
		while alseq[al_index] == '-':
			al_index += 1
		if seq_index == dr_list[curr_dr_num][0]:
			start_index = al_index
		if seq_index == dr_list[curr_dr_num][1] - 1:
			end_index = al_index + 1
			alseq_dr_crds.append([start_index, end_index])
			curr_dr_num += 1
		seq_index += 1
	return alseq_dr_crds


### COUNTING NUMBERS OF VARIOUS AMINOACIDS IN PROTEOME IN ORDERED/DISORDERED regions
def make_fastadict(fafile_name):
	fastadict = dict()
	with open(fafile_name) as inhandle:
		for record in SeqIO.parse(inhandle, "fasta"):
			s = ""
			for c in record.seq:
				s += c
			fastadict[record.id.split('|')[2]] = s
	return fastadict

def aa_num_calc(proteome_fasta, dr_filename, psite_filename, outfile_name):
	psite_dict = make_psite_dict(psite_filename)
	dr_dict = build_dr_dict(dr_filename)
	fastadict = make_fastadict(proteome_fasta)
	aa_num_dict = dict()
	for seq_id in fastadict.keys():
		seq = fastadict[seq_id]
		for dr_crd in dr_dict[seq_id]:
			for i in range(dr_crd[0], dr_crd[1]):
				aacid = seq[i]
				if psite_dict.get(seq_id) and psite_dict[seq_id].get(i):
					aacid = 'p' + aacid
				if aa_num_dict.get(aacid):
					aa_num_dict[aacid] += 1
				else:
					aa_num_dict[aacid] = 1
	with open(outfile_name, 'w') as outhandle:
		for k in aa_num_dict.keys():
			outhandle.write("{}\t{}\n".format(k, aa_num_dict[k]))



