"""
filter_sites

Filters the initial list of sites and returns the list without phosphosites
present only in protein isoforms or in proteins in which listed phosphosites
did not correspond to listed aminoacid. 

input:
1. Psite_table (reference)
2. Psite_table (being filtered)
3. Proteome_file

output:
Psite table with filtered phosphosites
"""

from Bio import SeqIO
from optparse import OptionParser


def make_psite_dict(psite_filename_ref):
	psite_num = 0
	psite_dict = dict()
	id_dict = dict()
	with open(psite_filename_ref) as psites:
		for psite_rec in psites:
			if psite_rec.startswith("UniProtID"):
				continue
			psite_rec = psite_rec.strip().split()
			if '-' in psite_rec[1]:
				continue
			psite_num += 1
			id_dict[psite_rec[0]] = psite_rec[1]
			if psite_dict.get(psite_rec[0]):
				psite_dict[psite_rec[0]][eval(psite_rec[2]) - 1] = psite_rec[3]
			else:
				psite_dict[psite_rec[0]] = dict()
				psite_dict[psite_rec[0]][eval(psite_rec[2]) - 1] = psite_rec[3]
#	calc = 0
#	for k in psite_dict.keys():
#		calc += 1
#	print calc
	return psite_dict, id_dict, psite_num


def find_incons_proteins(proteome_fasta, psite_dict, id_dict):
	inc_num = 0
	d = 0
	inc_prots = dict()
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
						inc_prots[seq_id] = True
#						print seq_id, crd, psite_dict[seq_id][crd], record.seq[crd]
				except IndexError:
					inc_num += 1
					d += 1
					inc_prots[seq_id] = True
	return inc_prots
#	print inc_num, d
#	print set(inc_prots)


def filter_psites(psite_dict, inc_prots, psite_filename_filt):
	with open(psite_filename_filt) as inhandle:
		for s in inhandle:
			str_list = s.strip().split()
			if str_list[0] == "UniProtID":
				print s.strip()
				continue
			if inc_prots.get(str_list[0]):
				continue
			if psite_dict.get(str_list[0]):
				if psite_dict[str_list[0]].get(eval(str_list[2]) - 1):
					print s.strip()


parser = OptionParser()
parser.add_option("-p", "--proteome_fasta", help="Proteome_file")
parser.add_option("-1", "--psite_filename_ref", help="Psite_table (reference)")
parser.add_option("-2", "--psite_filename_filt", help="Psite_table (being filtered)")
opt, args = parser.parse_args()


psite_dict, id_dict, psite_num = make_psite_dict(opt.psite_filename_ref)
#print psite_num
inc_prots = find_incons_proteins(opt.proteome_fasta, psite_dict, id_dict)
filter_psites(psite_dict, inc_prots, opt.psite_filename_filt)




