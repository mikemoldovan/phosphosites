"""
prepare_muttable_for_permuts

1. Leaves only mutations from alignment columns with a given aminoacid
2. Assigns phosphosites
3. Enimerates all mutations

Input:
1. Aminoacid
2. Mutation table
3. Phosphosite table
4. Organism ID
5. Ordered/disordered region flag
6. File with disordered regions
7. PAML output directory

Output:
Processed mutation table
"""

from optparse import OptionParser
from Bio import SeqIO


mut_counter = 0


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
	for c in seq:
		if c == '-':
			align_crd += 1
			continue
		if psite_dict[uniprot_id].get(seq_crd):
			psite_crds[align_crd] = True
		align_crd += 1
		seq_crd += 1
	return psite_crds


def leave_specific_cols(gene_substitutions, aminoacid):

	global mut_counter

	genesubst_dict = dict()
	for mut_record in gene_substitutions:
		if genesubst_dict.get(mut_record[0]):
			genesubst_dict[mut_record[0]].append(mut_record)
		else:
			genesubst_dict[mut_record[0]] = [mut_record]

	gene_substitutions_new = []
	mut_counter_list = []

	for k in genesubst_dict.keys():
		for mut_record in genesubst_dict[k]:
			if (mut_record[5] == aminoacid) or (mut_record[6] == aminoacid):
				for mut_record_2 in genesubst_dict[k]:
					gene_substitutions_new.append(mut_record_2)
					mut_counter_list.append(mut_counter)
				mut_counter += 1
				break

	return gene_substitutions_new, mut_counter_list


def block_analyze(gene_substitutions, 
				  paml_run_outdir, 
				  disord_reg_dict, 
				  psite_dict, 
				  uniprot_id, 
				  org_id,
				  flag,
				  aminoacid):
	
	psite_crds = align_mark_psites(psite_dict, org_id, uniprot_id, paml_run_outdir)
	block_paral_mat = dict()
	gene_substitutions, mut_counter_list = leave_specific_cols(gene_substitutions, aminoacid)

	for i in range(len(gene_substitutions)):

		mut_record = gene_substitutions[i]
		mut_counter = mut_counter_list[i]

		is_disord_reg_anc = False
		is_disord_reg_des = False

		c_anc = mut_record[5]
		c_des = mut_record[6]
		alpos = eval(mut_record[0])
		
		if psite_crds.get(eval(mut_record[0])):
			if c_anc == aminoacid:
				c_anc = 'p' + c_anc
				mut_record[5] = 'p' + mut_record[5]
			if c_des == aminoacid:
				c_des = 'p' + c_des
				mut_record[6] = 'p' + mut_record[6]

		if flag == "all":
			print str(mut_counter) + '\t' + '\t'.join(mut_record)
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
			print str(mut_counter) + '\t' + '\t'.join(mut_record)
		if flag == 'o' and not is_disord_reg_anc and not is_disord_reg_des:
			print str(mut_counter) + '\t' + '\t'.join(mut_record)


def main(muttable_filename, 
		 psite_filename, 
		 disord_reg_filename, 
		 org_id, 
		 flag,
		 aminoacid,
		 paml_run_outdir):
	psite_dict = make_psite_dict(psite_filename)

	if flag != "all":
		disord_reg_dict = make_disord_reg_dict(disord_reg_filename)
	else:
		disord_reg_dict = dict()

	gene_substitutions = []

	with open(muttable_filename) as muttable:
		for s in muttable:
			s = s.strip().split()
			if s[0] == "pos_in_al":
				print "pos_id" + '\t' + '\t'.join(s)
				continue
			if s[0].startswith('#'):
				if not gene_substitutions:
					uniprot_id = s[0][1:]
					print s[0]
				else:
					block_analyze(gene_substitutions,
								  paml_run_outdir, 
								  disord_reg_dict, 
								  psite_dict, 
								  uniprot_id, 
								  org_id,
								  flag,
								  aminoacid)
					gene_substitutions = []
				uniprot_id = s[0][1:]
				print s[0]
				continue
			gene_substitutions.append(s)
		block_analyze(gene_substitutions,
					  paml_run_outdir, 
					  disord_reg_dict, 
					  psite_dict, 
					  uniprot_id, 
					  org_id,
					  flag,
					  aminoacid)


parser = OptionParser()
parser.add_option("-m", "--muttable_filename", help="print_mutations output table")
parser.add_option("-p", "--psite_filename", help="Phosphosite table")
parser.add_option("-d", "--disord_reg_filename", help="File with disordered regions")
parser.add_option("-r", "--paml_run_outdir", help="PAML output directory")
parser.add_option("-i", "--org_id", help="spec_id for which phosphosites are counted")
parser.add_option("-f", "--flag", help="flag: 'all' for all phosphosites, 'd' for disordered regions, 'o' for ordered")
parser.add_option("-a", "--aminoacid", help="Aminoacid ('S', 'T' or 'Y')")
opt, args = parser.parse_args()


main(opt.muttable_filename, 
	 opt.psite_filename, 
	 opt.disord_reg_filename, 
	 opt.org_id, 
	 opt.flag,
	 opt.aminoacid,
	 opt.paml_run_outdir)

