"""
cat_phospholists

Concatenates several list of phosphosites

Input:
Lists of phosphosites (comma separated)

Output:
Concatenated list
"""

from optparse import OptionParser
from disord_reg_adlib import *

def main(psite_file_list):
	common_dict = dict()
	for psite_list in psite_file_list.split(','):
		psite_dict = make_psite_dict2(psite_list)
		for uniprot_id in psite_dict.keys():
			if common_dict.get(uniprot_id):
				for crd in psite_dict[uniprot_id].keys():
					if not common_dict[uniprot_id].get(crd):
						common_dict[uniprot_id][crd] = psite_dict[uniprot_id][crd]
			else:
				common_dict[uniprot_id] = psite_dict[uniprot_id]

	for uniprot_id in common_dict.keys():
		for crd in common_dict[uniprot_id].keys():
			print '\t'.join(common_dict[uniprot_id][crd][:5])


parser = OptionParser()
parser.add_option("-i", "--psite_file_list", help="Lists of phosphosites (comma separated)")
opt, args = parser.parse_args()

main(opt.psite_file_list)