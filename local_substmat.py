"""
local_substmat

Builds two substitution matrices: the one around phosphosites and the one around
non-phosphorilated aminoacids

Input:
1. List of phosphosites
2. List of phosphorilated aminoacids
3. Radius: half the window size for the substitution matrix calculation
4. Mutation list (print_mutations output)
5. Proteome FASTA
6. Flag: in disordered region/not in disordered region yes/no/all (default "yes")
"""

from optparse import OptionParser
from disord_reg_adlib import *


def mutdict_append(mutdict, phospho_aas, anc_let, des_let, psites, subst_pos):
	if psites.get(subst_pos):
		if anc_let in phospho_aas:
			anc_let = 'p' + anc_let
		if des_let in phospho_aas: 
			des_let = 'p' + des_let
	if mutdict.get(anc_let):
		if mutdict[anc_let].get(des_let):
			mutdict[anc_let][des_let] += 1
		else:
			mutdict[anc_let][des_let] = 1
	else:
		mutdict[anc_let] = dict()
		mutdict[anc_let][des_let] = 1


def main(mutlist_file, psite_filename, in_dr, phospho_aas, proteome_fasta, radius):
	psite_dict = make_psite_dict(psite_filename)
	proteome_dict = make_fastadict(proteome_fasta)
	mutdict_np = dict()
	mutdict_p = dict()
	with open(mutlist_file) as mutlist:
		for mutrec in mutlist:
			mutrec = mutrec.strip().split()
			if mutrec[0].startswith('#'):
				if len(mutrec) == 1:
					gene_id = mutrec[0][1:]
					if proteome_dict.get(gene_id):
						seq = proteome_dict[gene_id]
						seqlen = len(seq)
					else:
						continue
					if psite_dict.get(gene_id):
						psites = psite_dict[gene_id]
					else:
						psites = dict()
				continue
			if mutrec[8] == 'N':
				continue
			if mutrec[8] == -1:
				continue
			if in_dr == 'y' and mutrec[7] == '*':
				continue
			if in_dr == 'n' and mutrec[7] != '*':
				continue
			anc_let = mutrec[5]
			des_let = mutrec[6]
			if anc_let in phospho_aas or des_let in phospho_aas:
				continue
			is_in_rad = False
			is_phos_reg = False
			subst_pos = eval(mutrec[8])
			for i in (subst_pos - radius, subst_pos + radius + 1):
				if i < 0 or i >= seqlen:
					continue
				if seq[i] in phospho_aas:
					is_in_rad = True
					if psites.get(i):
						is_phos_reg = True
			if not is_in_rad:
				continue
			if is_phos_reg:
				mutdict_append(mutdict_p, phospho_aas, anc_let, des_let, psites, subst_pos)
			else:
				mutdict_append(mutdict_np, phospho_aas, anc_let, des_let, psites, subst_pos)
	print {"mutdict_p": mutdict_p, "mutdict_np": mutdict_np}


parser = OptionParser()
parser.add_option("-m", "--mutlist_file", help="Mutation list (print_mutations output)")
parser.add_option("-p", "--psite_filename", help="Phosphosite table")
parser.add_option("-d", "--proteome_fasta", help="Proteome file (optonal, if there is no table with site numbers)", default = "_")
parser.add_option("-f", "--in_dr", help="Flag: in disordered region/not in disordered region y/n/all (default 'y')", default = 'y')
parser.add_option("-r", "--radius", help="Radius: half the window size for the substitution matrix calculation", default = 'y')
parser.add_option("-a", "--phospho_aas", help="List of phosphorilated aminoacids (default 'STY')", default = "STY")
opt, args = parser.parse_args()


main(opt.mutlist_file, 
	 opt.psite_filename, 
	 opt.in_dr, 
	 opt.phospho_aas, 
	 opt.proteome_fasta, 
	 eval(opt.radius))

