"""
substvec_permut_count

Count substitution vector for phosphorilated aminoacid and permutation vectors
for its unphosphorilated analogue. Works with trimmed mutation table
(prepare_muttable_for_permuts.py output)

Input:
1. Mutation table
2. Number of permutations
3. Job ID
4. Aminoacid
5. BLOSUM matrix

Output:
1. Table with permutations for (p)N1 -> N2 mutation vectors
2. Table with permutations for N1 -> (p)N2 mutation vectors 
3. Table with permutations for (p)N1 -> N2 H values 
"""


aacids = list("ACDEFGHIKLMNPQRSTVWY")


import random
from optparse import OptionParser
from os import remove


def make_mut_dict(muttable_filename):
	curr_pos_id = -1
	pos_list = []
	p_pos_list = []
	mut_dict = dict()

	mutlist = []
	ispsite = False

	with open(muttable_filename) as muttable:
		for s in muttable:
			if s.startswith('#'):
				curr_pos_id = -1
				continue
			if s.startswith("pos_id"):
				continue
			s = s.strip().split()
			pos_id = eval(s[0])
			if pos_id != curr_pos_id and mutlist:
				if ispsite:
					p_pos_list.append(curr_pos_id)
				else:
					pos_list.append(curr_pos_id)
				mut_dict[curr_pos_id] = mutlist
				mutlist = []
				curr_pos_id = pos_id
				ispsite = False
			mutlist.append((s[6],s[7]))
			if s[6][0] == 'p' or s[7][0] == 'p':
				ispsite = True
		if mutlist:
			if ispsite:
				p_pos_list.append(curr_pos_id)
			else:
				pos_list.append(curr_pos_id)
			mut_dict[curr_pos_id] = mutlist

	return mut_dict, pos_list, p_pos_list


def infer_blosum_mat(blosum_mat_filename):
	blosum_mat = dict()
	with open(blosum_mat_filename) as bl_mat:
		for s in bl_mat:
			s = s.strip().split()
			if s[0] == '#':
				l = len(s)
				let_list = s
				for c in s[1:]:
					blosum_mat[c] = dict()
				continue
			for i in range(1, l):
				blosum_mat[s[0]][let_list[i]] = eval(s[i])
	return blosum_mat


def substvec_count(mut_dict,
				   p_pos_list,
				   from_psite_mutfile,
				   to_psite_mutfile,
				   from_psite_paralfile,
				   permut_id,
				   blosum_mat,
				   aacid):
	
	global aacids

	from_psite_substvec = dict()
	to_psite_substvec = dict()
	from_psite_paralvec = dict()

	for c in aacids:
		from_psite_substvec[c] = 0
		to_psite_substvec[c] = 0
		from_psite_paralvec[c] = 0

	for pos_id in p_pos_list:
		from_mutlist = []
		for mutation in mut_dict[pos_id]:
			if mutation[0][-1] == aacid:
				from_psite_substvec[mutation[1]] += 1
				from_mutlist.append(mutation)
			if mutation[1][-1] == aacid:
				to_psite_substvec[mutation[0]] += 1
		l = len(from_mutlist)
		for i in range(l):
			for j in range(i + 1, l):
				from_psite_paralvec[from_mutlist[i][1]] += blosum_mat[from_mutlist[i][1][-1]][from_mutlist[j][1][-1]]

	from_psite_str_vec = '\t'.join([str(from_psite_substvec[c]) for c in aacids])
	to_psite_str_vec = '\t'.join([str(to_psite_substvec[c]) for c in aacids])
	from_psite_str_paralvec = '\t'.join([str(from_psite_paralvec[c]) for c in aacids])

	from_psite_mutfile.write(str(permut_id) + '\t' + from_psite_str_vec + '\n')
	to_psite_mutfile.write(str(permut_id) + '\t' + to_psite_str_vec + '\n')
	from_psite_paralfile.write(str(permut_id) + '\t' + from_psite_str_paralvec + '\n')


def mutcount(mut_dict, p_pos_list, flag, aacid):
	mutnum = 0
	for pos_id in p_pos_list:
		for mut in mut_dict[pos_id]:
			if flag == "from":
				if mut[0][-1] == aacid:
					mutnum += 1
			if flag == "to":
				if mut[1][-1] == aacid:
					mutnum += 1
	return mutnum


def main(muttable_filename, permut_num, job_id, aacid, blosum_mat_filename, flag):

	global aacids

	aacids = list(aacids)
	aacids.append('p' + aacid)

	mut_dict, pos_list, p_pos_list = make_mut_dict(muttable_filename)
	blosum_mat = infer_blosum_mat(blosum_mat_filename)

	from_psite_mutfile = open("{}from_psite_muts.txt".format(job_id), 'w')
	to_psite_mutfile = open("{}to_psite_muts.txt".format(job_id), 'w')
	from_psite_paralfile = open("{}from_psite_paral.txt".format(job_id), 'w')

	from_psite_mutfile.write('#' + '\t' + '\t'.join(aacids) + '\n')
	to_psite_mutfile.write('#' + '\t' + '\t'.join(aacids) + '\n')
	from_psite_paralfile.write('#' + '\t' + '\t'.join(aacids) + '\n')

	permut_id = -1

	substvec_count(mut_dict, 
				   p_pos_list, 
				   from_psite_mutfile, 
				   to_psite_mutfile, 
				   from_psite_paralfile,
				   permut_id,
				   blosum_mat,
				   aacid)

	psite_num = len(p_pos_list)
	psite_mutnum = mutcount(mut_dict, p_pos_list, flag, aacid)
#	print psite_mutnum
	permut_mutnums = []

	for i in range(100):
		permut_p_pos_list = random.sample(pos_list, psite_num)
		permut_mutnums.append(mutcount(mut_dict, permut_p_pos_list, flag, aacid))

	mean_permut_mutnum = sum(permut_mutnums)/100
#	print mean_permut_mutnum
	randsite_num = (psite_num*psite_mutnum)/mean_permut_mutnum

	print "Randsite number:\t{}\t{}".format(randsite_num, float(psite_mutnum)/mean_permut_mutnum)

	for permut_id in range(permut_num):
		permut_p_pos_list = random.sample(pos_list, randsite_num)
		substvec_count(mut_dict, 
					   permut_p_pos_list, 
					   from_psite_mutfile, 
					   to_psite_mutfile, 
					   from_psite_paralfile,
					   permut_id,
					   blosum_mat,
					   aacid)

	from_psite_mutfile.close()
	to_psite_mutfile.close()
	from_psite_paralfile.close()

	if flag == "from":
		print "{}from_psite_muts.txt".format(job_id)
		print "{}from_psite_paral.txt".format(job_id)
		remove("{}to_psite_muts.txt".format(job_id))
	else:
		print "{}to_psite_muts.txt".format(job_id)
		remove("{}from_psite_muts.txt".format(job_id))
		remove("{}from_psite_paral.txt".format(job_id))


parser = OptionParser()
parser.add_option("-m", "--muttable_filename", help="prepare_muttable_for_permuts output")
parser.add_option("-n", "--permut_num", help="Permutation number")
parser.add_option("-i", "--job_id", help="Job ID (default '')", default="")
parser.add_option("-a", "--aacid", help="Aminoacid ('S', 'T' or 'Y')")
parser.add_option("-b", "--blosum_mat_filename", help="BLOSUM matrix")
parser.add_option("-f", "--flag", help="Flag: 'from' or 'to'")
opt, args = parser.parse_args()


main(opt.muttable_filename, eval(opt.permut_num), opt.job_id, opt.aacid, opt.blosum_mat_filename, opt.flag)