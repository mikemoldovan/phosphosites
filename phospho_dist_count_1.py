"""
phospho_dist_count_1

yields the distribution of distances between phosphosites

Input:
1. Phosphorilated aminoacids (like 'S' or "YST")
2. FASTA file with protein sequences (proteome)
3. Phosphosite table
4. Number of permutations for the null distribution

Output:
1. Null distribution of distances
2. Real distribution of distances
"""


import random
from Bio import SeqIO
from optparse import OptionParser


def make_psite_dict(psite_filename, psite_aas):
	psite_dict = dict()
	with open(psite_filename) as psite_file:
		for s in psite_file:
			s = s.strip().split()
			if s[0] == "UniProtID":
				continue
			if psite_dict.get(s[0]):
				if s[3] in psite_aas:
					psite_dict[s[0]][eval(s[2])] = s[3]
			else:
				if s[3] in psite_aas:
					psite_dict[s[0]] = dict()
					psite_dict[s[0]][eval(s[2])] = s[3]
	return psite_dict


def make_null_dist(seq, psite_aas, permut_num, psite_num, null_outfile):
	psite_aas_poss = []
	for i in range(len(seq)):
		if seq[i] in psite_aas:
			psite_aas_poss.append(i)
	for i in range(permut_num):
		poss = sorted(random.sample(psite_aas_poss, psite_num))
		for i in range(1, psite_num):
			null_outfile.write(str(poss[i] - poss[i - 1]) + '\n')


def main(psite_aas, protein_fasta, psite_filename, permut_num):
	psite_aas = list(psite_aas)
	null_outfile = open("null_dist.txt", 'w')
	outfile = open("psite_distance_dist.txt", 'w')
	psite_dict = make_psite_dict(psite_filename, psite_aas)
	with open(protein_fasta) as proteome_fasta:
		for record in SeqIO.parse(proteome_fasta, "fasta"):
			_id = record.id.split('|')[2]
			if not psite_dict.get(_id):
				continue
			seq = ""
			for c in record.seq:
				seq += c
			psite_crds = list(psite_dict[_id].keys())
			psite_crds = sorted(psite_crds)
			psite_crds_new = []
			for crd in psite_crds:
				try:
					if seq[crd - 1] in psite_aas:
						psite_crds_new.append(crd)
				except:
					break
			len_psite_crds = len(psite_crds_new)
			if len_psite_crds < 2:
				continue
			try:
				make_null_dist(seq, psite_aas, permut_num, len_psite_crds, null_outfile)
			except:
				print psite_crds_new, seq, _id
			for i in range(1, len_psite_crds):
				outfile.write(str(psite_crds_new[i] - psite_crds_new[i - 1]) + '\n')

	null_outfile.close()
	outfile.close()
	print "null distribution in: null_dist.txt"
	print "psite distance distribution in: psite_distance_dist.txt"


parser = OptionParser()
parser.add_option("-a", "--psite_aas", help="Phosphorilated aminoacids (like 'S' or 'YST')")
parser.add_option("-p", "--protein_fasta", help="FASTA file with protein sequences (proteome)")
parser.add_option("-s", "--psite_filename", help="Phosphosite table")
parser.add_option("-n", "--permut_num", help="Number of permutations for the null distribution", default='2')
opt, args = parser.parse_args()


main(opt.psite_aas, opt.protein_fasta, opt.psite_filename, eval(opt.permut_num))