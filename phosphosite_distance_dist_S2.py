"""
phosphosite_distance_dist_S2

Count distributions of distances
between phosphosites for different disordered regions

Input:
1. Table with phosphosites
2. Table with disordered regions
3. File with proteome
4. Number of permutations 
5. Aminoacids for which phosphosites to be counted (default 'STY')

Output:

1. Table with distribution of distances between phosphosites
2. Table with permutations for [2]
"""


import random
from optparse import OptionParser
#from Bio import SeqIO
from disord_reg_adlib import *


def calc_s2(dr_dict, psite_dict, outfile):
	num_dr_sites = dict()
	for gene_id in dr_dict.keys():
		if not psite_dict.get(gene_id):
			continue
		psite_crds = sorted(list(psite_dict[gene_id].keys()))
#		print psite_crds
		for crd_pair in dr_dict[gene_id]:
#			print crd_pair
			in_crd_pair = []
			for psite_crd in psite_crds:
				if psite_crd >= crd_pair[0] and psite_crd < crd_pair[1]:
					phos_let = psite_dict[gene_id][psite_crd]
					if num_dr_sites.get(phos_let):
						num_dr_sites[phos_let] += 1
					else:
						num_dr_sites[phos_let] = 1
					in_crd_pair.append(psite_crd)
			if not in_crd_pair:
				continue
			for i in range(1, len(in_crd_pair)):
				crd_str = str(crd_pair[0]) + '-' + str(crd_pair[1])
				s2 = in_crd_pair[i] - in_crd_pair[i-1]
#				outfile.write("{}\t{}\t{}\t{}\t{}\n".format(gene_id,
#											crd_str,
#											in_crd_pair[i-1],
#											in_crd_pair[i],
#											s2))
				outfile.write(str(s2) + '\n')
	print(num_dr_sites)
	return num_dr_sites


def obtain_rand_s2(fastadict, dr_dict, num_dr_sites, phospho_aas):
	site_ids = dict()
	all_sitedict = dict()
	for aa in phospho_aas:
		all_sitedict[aa] = dict()
	for aa in phospho_aas:
		site_ids[aa] = 0
	for gene_id in fastadict.keys():
		if not dr_dict.get(gene_id):
			continue
		for dr in dr_dict[gene_id]:
			for crd in range(dr[0],dr[1]):
				c = fastadict[gene_id][crd]
				if c in phospho_aas:
					all_sitedict[c][site_ids[c]] = "{}|{}".format(gene_id, crd)
					site_ids[c] += 1
	print(site_ids)
	rand_psite_dict = dict()
	for c in phospho_aas:
#		print all_sitedict[c]
		for sitenum in range(num_dr_sites[c]):
			diff_val = False
			while diff_val == False:
				randval = random.randint(0, site_ids[c] - 1)
				psite = all_sitedict[c][randval].split('|')
				gene_id = psite[0]
				crd = eval(psite[1])
				if rand_psite_dict.get(gene_id):
					if not rand_psite_dict[gene_id].get(crd):
						rand_psite_dict[gene_id][crd] = c
						diff_val = True
				else:
					rand_psite_dict[gene_id] = dict()
					rand_psite_dict[gene_id][crd] = c
					diff_val = True
	return rand_psite_dict


def main(psite_filename,
		 dr_filename,
		 proteome_fasta,
		 permut_num,
		 phospho_aas,
		 job_id):

	psite_dict = make_psite_dict(psite_filename)
	dr_dict = build_dr_dict(dr_filename)
	fastadict = make_fastadict(proteome_fasta)

	S2_outfile = open("{}S2_stats.txt".format(job_id), 'w')
	S2_permuts = open("{}S2_permuts.txt".format(job_id), 'w')

	num_dr_sites = calc_s2(dr_dict, psite_dict, S2_outfile)
	for i in range(permut_num):
		rand_psite_dict = obtain_rand_s2(fastadict, dr_dict, num_dr_sites, phospho_aas)
		calc_s2(dr_dict, rand_psite_dict, S2_permuts)
	S2_outfile.close()
	S2_permuts.close()


parser = OptionParser()
parser.add_option("-p", "--psite_filename", help="Table with phosphosites")
parser.add_option("-r", "--dr_filename", help="Table with disordered regions")
parser.add_option("-f", "--proteome_fasta", help="File with proteome")
parser.add_option("-n", "--permut_num", help="Number of permutations (default='2')", default='2')
parser.add_option("-a", "--phospho_aas", help="Aminoacids for which phosphosites to be counted (default 'STY')", default="STY")
parser.add_option("-j", "--job_id", help="Job identifier (default '')", default='')
opt, args = parser.parse_args()


main(opt.psite_filename,
	 opt.dr_filename,
	 opt.proteome_fasta,
	 eval(opt.permut_num),
	 opt.phospho_aas,
	 opt.job_id)

