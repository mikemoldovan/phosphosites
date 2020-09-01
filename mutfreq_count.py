"""
mutfreq_count

Counts mutational frequencies versus phosphorilated and non-phosphorilated aminoacids.
Performs Fisher's exact test to check for significance

Input:
1. Mutation list (print_mutations output)
2. Phosphosite table
3. Proteome file (optonal, if there is no table with site numbers)
4. Disordered region file (optonal, if there is no table with site numbers)
5. Table with site numbers (optional)
6. Flag: in disordered region/not in disordered region yes/no/all (default "yes")
7. List of phosphorilated aminoacids (default "STY")
"""

import pickle
from scipy import stats
from optparse import OptionParser
from disord_reg_adlib import *
from sys import stdout


def sitenum_file_parse(sitenum_filename, phospho_aas):
	sitenum_dict = dict()
	with open(sitenum_filename) as sitenums:
		for s in sitenums:
			s = s.strip().split()
			if s[0][0] == 'p' and s[0][-1] in phospho_aas:
				sitenum_dict[s[0]] = eval(s[1])
			elif s[0][0] == 'p' and s[0][-1] not in phospho_aas:
				sitenum_dict[s[0][-1]] = eval(s[1])
			else:
				sitenum_dict[s[0]] = eval(s[1])
	return sitenum_dict

"""
def print_results(mutdict, sitenum_filename, phospho_aas):
	aa_list = list("ACEDGFIHKMLNQPSRTWVY")
	for aacid in phospho_aas:
		aa_list.append('p' + aacid)
	for aa in phospho_aas:
		print aa
		p_aa = 'p' + aa
		for aa_des in aa_list:
			if mutdict[aa].get(aa_des):
				ref_mutnum = mutdict[aa][aa_des]
			else:
				ref_mutnum = 0
			all_ref_muts = mutdict_colsum(mutdict, aa)
			if mutdict[p_aa].get(aa_des):
				mutnum = mutdict[p_aa][aa_des]
			else:
				mutnum = 0
			all_muts = mutdict_colsum(mutdict, p_aa)
			oddsratio, pvalue = stats.fisher_exact([[mutnum, all_muts], [ref_mutnum, all_ref_muts]])
			print aa, aa_des, mutnum, all_muts, ref_mutnum, all_ref_muts, float(mutnum)/all_muts, float(ref_mutnum)/all_ref_muts, pvalue
"""

def main(mutlist_file, psite_filename, proteome_fasta, dr_filename, sitenum_filename, in_dr, phospho_aas):
	if sitenum_filename == 'n':
		sitenum_filename = "sitenums.txt"
		aa_num_calc(proteome_fasta, dr_filename, psite_filename, sitenum_filename)
	mutdict = dict()
	psite_dict = make_psite_dict(psite_filename)
	with open(mutlist_file) as mutlist:
		for mutrec in mutlist:
			mutrec = mutrec.strip().split()
			if mutrec[0].startswith('#'):
				if len(mutrec) == 1:
					gene_id = mutrec[0][1:]
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
			if psites.get(eval(mutrec[8])):
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
	sitenum_dict = sitenum_file_parse(sitenum_filename, phospho_aas)
	outdict = {"mutdict": mutdict, "sitenumdict": sitenum_dict}
	print(outdict)


parser = OptionParser()
parser.add_option("-m", "--mutlist_file", help="Mutation list (print_mutations output)")
parser.add_option("-p", "--psite_filename", help="Phosphosite table")
parser.add_option("-d", "--proteome_fasta", help="Proteome file (optonal, if there is no table with site numbers)", default = "_")
parser.add_option("-r", "--dr_filename", help="Disordered region file (optonal, if there is no table with site numbers)", default = "_")
parser.add_option("-i", "--sitenum_filename", help="Table with site numbers (optional), 'n' for no file", default = "_")
parser.add_option("-f", "--in_dr", help="Flag: in disordered region/not in disordered region y/n/all (default 'y')", default = 'y')
parser.add_option("-a", "--phospho_aas", help="List of phosphorilated aminoacids (default 'STY')", default = "STY")
opt, args = parser.parse_args()

main(opt.mutlist_file, 
	 opt.psite_filename, 
	 opt.proteome_fasta, 
	 opt.dr_filename, 
	 opt.sitenum_filename, 
	 opt.in_dr, 
	 opt.phospho_aas)


"""
python ../scripts/mutfreq_count.py -m human_oggs_mammalia_paml_mutations2.txt -p sitelists_filtered/human_iptmnet.tsv -d proteomes/human_proteome.fa -r human_proteome_disordered_regs.txt -i n > ../mutfreq_count/human_np_iptmnet.txt
python ../scripts/mutfreq_count.py -m human_oggs_mammalia_paml_mutations2.txt -p sitelists_filtered/human_iptmnet_minscore2.tsv -d proteomes/human_proteome.fa -r human_proteome_disordered_regs.txt -i n > ../mutfreq_count/human_np_iptmnet_minscore2.txt
python ../scripts/mutfreq_count.py -m human_oggs_mammalia_paml_mutations2.txt -p sitelists_filtered/human_iptmnet_minscore3.tsv -d proteomes/human_proteome.fa -r human_proteome_disordered_regs.txt -i n > ../mutfreq_count/human_np_iptmnet_minscore3.txt
python ../scripts/mutfreq_count.py -m human_oggs_mammalia_paml_mutations2.txt -p sitelists_filtered/human_mouse_cons_iptmnet.tsv -d proteomes/human_proteome.fa -r human_proteome_disordered_regs.txt -i n > ../mutfreq_count/human_np_mouse_cons.txt
"""