"""
uniprot_download

download protein sequences for a given phosphosite file

input:
phosphosite_file

output:
FASTA sequences of phosphorilated proteins
"""

import urllib
from optparse import OptionParser

def make_id_dict(phosphosite_file):
	id_dict = dict()
	with open(phosphosite_file) as pfile:
		for s in pfile:
			s = s.strip().split()
			if s[0] == "UniProtID":
				continue
			id_dict[s[0]] = True
	return id_dict


def download_proteins(id_dict):
	for _id in id_dict.keys():
		link = "https://www.uniprot.org/uniprot/{}.fasta".format(_id)
		f = urllib.urlopen(link)
		myfile = f.read()
		print myfile


def main(phosphosite_file):
	id_dict = make_id_dict(phosphosite_file)
	download_proteins(id_dict)


parser = OptionParser()
parser.add_option("-i", "--psite_info_filename", help="file in format: niprotID\IPTMnetID\position\ residue\score (iptmnet_info output)")
opt, args = parser.parse_args()

main(opt.psite_info_filename)