"""
substmat_count_2

Count substitution matrix from the print_mutations output table

Input:
1. print_mutations output table
2. Phosphosite table
3. File with disordered regions
4. spec_id
5. PAML output directory

Output:
substitution matrix
"""

from optparse import OptionParser
from Bio import SeqIO


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


def make_disord_reg_dict(disord_reg_filename):
	disord_reg_dict = dict()
	with open(disord_reg_filename) as disord_regs:
		for protein_rec in disord_regs:
			protein_rec = protein_rec.strip().split()
			try:
				a = protein_rec[1]
			except:
				print protein_rec
			if disord_reg_dict.get(protein_rec[0]):
				disord_reg_dict[protein_rec[0]][protein_rec[1]] = []
			else:
				disord_reg_dict[protein_rec[0]] = dict()
				disord_reg_dict[protein_rec[0]][protein_rec[1]] = []
			if len(protein_rec) < 3:
				disord_reg_dict[protein_rec[0]][protein_rec[1]] = []
				continue
			for disord_reg in protein_rec[2].split('|'):
				reg = disord_reg.split('-')
				reg = [eval(reg[0]) - 1, eval(reg[1])]
				disord_reg_dict[protein_rec[0]][protein_rec[1]].append(reg)
	return disord_reg_dict


def obtain_seq(faflie_name, org_id):
	with open(faflie_name) as infasta:
		for record in SeqIO.parse(infasta, "fasta"):
			if record.id.split('_')[-1] == org_id:
				seq = ""
				for c in record.seq:
					seq += c
				return seq


def align_mark_psites(psite_dict, org_id, uniprot_id, paml_run_outdir):
	seq = obtain_seq(paml_run_outdir + uniprot_id + ".fa", org_id)
	psite_crds = dict()
	align_crd = 0
	seq_crd = 0
	if not seq:
		return psite_crds
	if not psite_dict.get(uniprot_id):
		return psite_crds
#	print "alseq", seq
	for c in seq:
		if c == '-':
			align_crd += 1
			continue
		if psite_dict[uniprot_id].get(seq_crd):
			psite_crds[align_crd] = True
		align_crd += 1
		seq_crd += 1
#	print "psite_crds", psite_crds
	return psite_crds


def substmat_append(substmat, c_anc, c_des):
	if substmat.get(c_anc):
		if substmat[c_anc].get(c_des):
			substmat[c_anc][c_des] += 1
		else:
			substmat[c_anc][c_des] = 1
	else:
		substmat[c_anc] = dict()
		substmat[c_anc][c_des] = 1

	
def block_analyze_disord(substmat,
						 gene_substitutions, 
						 paml_run_outdir, 
						 disord_reg_dict, 
						 psite_dict, 
						 uniprot_id, 
						 org_id,
						 flag):
#	print "\n\n"
	psite_crds = align_mark_psites(psite_dict, org_id, uniprot_id, paml_run_outdir)
#	try:
#		print "uniprot_id", uniprot_id
#		print "psites_original", psite_dict[uniprot_id]
#		print "psite_crds", psite_crds
#		print "disord_regs", disord_reg_dict[uniprot_id]
#	except:
#		pass
	for mut_record in gene_substitutions:
		is_disord_reg_anc = False
		is_disord_reg_des = False

		c_anc = mut_record[5]
		c_des = mut_record[6]
		
		if psite_crds.get(eval(mut_record[0])):
			if c_anc in ('S','T','Y'):
				c_anc = 'p' + c_anc
			if c_des in ('S','T','Y'):
				c_des = 'p' + c_des

		if flag == "all":
			substmat_append(substmat, c_anc, c_des)
			continue

		try:
			disord_regs_anc = disord_reg_dict[uniprot_id][mut_record[1]]
			disord_regs_des = disord_reg_dict[uniprot_id][mut_record[2]]
		except:
			if flag == 'o':
				disord_regs_anc = [[0,0]]
				disord_regs_des = [[0,0]]
			else:
				continue
		
		pos_anc = eval(mut_record[3])
		pos_des = eval(mut_record[4])
		
		for crd_pair in disord_regs_anc:
			if pos_anc >= crd_pair[0] and pos_anc < crd_pair[1]:
				is_disord_reg_anc = True
				break

		for crd_pair in disord_regs_des:
			if pos_des >= crd_pair[0] and pos_des < crd_pair[1]:
				is_disord_reg_des = True
				break

		if flag == 'd' and is_disord_reg_anc and is_disord_reg_des:
			substmat_append(substmat, c_anc, c_des)
		if flag == 'o' and not is_disord_reg_anc and not is_disord_reg_des:
			substmat_append(substmat, c_anc, c_des)


def print_subst_mat(subst_mat):
	colnames = sorted(list(subst_mat.keys()))
	colnames_l = len(colnames)
	rownames = []
	for k1 in subst_mat.keys():
		for k2 in subst_mat[k1].keys():
			rownames.append(k2)
	rownames = sorted(list(set(rownames)))
	rownames_l = len(rownames)
	print '\t'.join(['#'] + colnames)
	for c1 in colnames:
		row = [c1]
		for c2 in rownames:
			if subst_mat[c1].get(c2):
				row.append(subst_mat[c1][c2])
			else:
				row.append(0)
		print '\t'.join(map(str, row))


def main(muttable_filename, 
		 psite_filename, 
		 disord_reg_filename, 
		 paml_run_outdir, 
		 org_id, 
		 flag):
	psite_dict = make_psite_dict(psite_filename)
#	print "psite_dict built"
	if flag != "all":
		disord_reg_dict = make_disord_reg_dict(disord_reg_filename)
	else:
		disord_reg_dict = dict()
#	print "disord_reg_dict built, counting substmat"
	substmat = dict()
	gene_substitutions = []

	with open(muttable_filename) as muttable:
		for s in muttable:
			s = s.strip().split()
			if s[0] == "pos_in_al":
				continue
			if s[0].startswith('#'):
				if not gene_substitutions:
					uniprot_id = s[0][1:]
				else:
					block_analyze_disord(substmat,
										 gene_substitutions,
										 paml_run_outdir, 
										 disord_reg_dict, 
										 psite_dict, 
										 uniprot_id, 
										 org_id,
										 flag)
#					print gene_substitutions
					gene_substitutions = []
#					print "psites", psite_dict[uniprot_id]
#					print "uniprot_id", uniprot_id
#					break
				uniprot_id = s[0][1:]
				continue
			gene_substitutions.append(s)
		block_analyze_disord(substmat,
							 gene_substitutions,
							 paml_run_outdir, 
							 disord_reg_dict, 
							 psite_dict, 
							 uniprot_id, 
							 org_id,
							 flag)
	print_subst_mat(substmat)
	print substmat


parser = OptionParser()
parser.add_option("-m", "--muttable_filename", help="print_mutations output table")
parser.add_option("-p", "--psite_filename", help="Phosphosite table")
parser.add_option("-d", "--disord_reg_filename", help="File with disordered regions")
parser.add_option("-r", "--paml_run_outdir", help="PAML output directory")
parser.add_option("-i", "--org_id", help="spec_id for which phosphosites are counted")
parser.add_option("-f", "--flag", help="flag: 'all' for all phosphosites, 'd' for disordered regions, 'o' for ordered")
opt, args = parser.parse_args()


main(opt.muttable_filename, 
	 opt.psite_filename, 
	 opt.disord_reg_filename, 
	 opt.paml_run_outdir, 
	 opt.org_id, 
	 opt.flag)

