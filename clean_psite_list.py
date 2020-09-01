"""
clean_orig_psite_list

input:
1. Psite_table
2. Proteome_file
"""

from Bio import SeqIO
from optparse import OptionParser


def make_psite_dict(psite_filename):
	psite_num = 0
	psite_dict = dict()
	id_dict = dict()
	with open(psite_filename) as psites:
		for psite_rec in psites:
			if psite_rec.startswith("UniProtID"):
				continue
			psite_rec = psite_rec.strip().split()
#			if '-' in psite_rec[1]:
#				continue
			psite_num += 1
			id_dict[psite_rec[0]] = psite_rec[1]
			if psite_dict.get(psite_rec[0]):
				psite_dict[psite_rec[0]][eval(psite_rec[2]) - 1] = psite_rec[3]
			else:
				psite_dict[psite_rec[0]] = dict()
				psite_dict[psite_rec[0]][eval(psite_rec[2]) - 1] = psite_rec[3]
	calc = 0
	for k in psite_dict.keys():
		calc += 1
	print "Total phosphoproteins: ", calc
	return psite_dict, id_dict, psite_num


def find_inconsistencies(proteome_fasta, psite_dict, id_dict):
	inc_num = 0
	d = 0
	inc_prots = []
	with open(proteome_fasta) as infasta:
		for record in SeqIO.parse(infasta, "fasta"):
			seq_id = record.id.split('|')[2]
			seq_id_2 = record.id.split('|')[1]
			try:
				if id_dict[seq_id] != seq_id_2:
					continue
			except:
				continue
			if not psite_dict.get(seq_id):
				continue
			for crd in psite_dict[seq_id].keys():
				try:
					a = record.seq[crd]
					if record.seq[crd] != psite_dict[seq_id][crd]:
						inc_num += 1
						inc_prots.append(seq_id)
#						print seq_id, crd, psite_dict[seq_id][crd], record.seq[crd]
				except IndexError:
					inc_num += 1
					d += 1
	print "Number of inconsistent proteins: ", inc_num
	print "IndexErrors: ", d
#	print set(inc_prots)


parser = OptionParser()
parser.add_option("-p", "--proteome_fasta", help="Proteome_file")
parser.add_option("-s", "--psite_filename", help="Psite_table")
opt, args = parser.parse_args()

psite_dict, id_dict, psite_num = make_psite_dict(opt.psite_filename)
print "Psite number: ", psite_num
find_inconsistencies(opt.proteome_fasta, psite_dict, id_dict)
