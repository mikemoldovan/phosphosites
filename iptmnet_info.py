"""
iptmnet_info

obtain information about all phosphorilation sites from iptmnet

input:
1. Tax_ID
2. Size of the single sample (default 100)
3. Maximal number of proteins

output:
File in format:
UniprotID\IPTMnetID\position\residue\score
"""

from os import system, remove, listdir
from optparse import OptionParser

def fileeval(filename):
	with open(filename) as eval_file:
		file_str = ""
		for string in eval_file:
			file_str += string.strip()
		file_str = file_str.replace("true", "True")
		file_str = file_str.replace("false", "False")
		file_str = file_str.replace("null", "None")
	return eval(file_str)


def outfile_write_strs(outfile, uniprot_id, iptmnet_id, ptm_dicts_dict):
	for k in ptm_dicts_dict.keys():
		ptm_dicts = ptm_dicts_dict[k]
		for ptm_dict in ptm_dicts:
			if ptm_dict["ptm_type"] != "Phosphorylation":
				continue
			outfile.write("{}\t{}\t{}\t{}\t{}\n".format(uniprot_id,
														k,
														ptm_dict["site"][1:],
														ptm_dict["residue"],
														ptm_dict["score"]))


def main(tax_id, sample_size, ptm_prot_num):
	for i in range(0, ptm_prot_num, sample_size):
		endpos = i + sample_size if i + sample_size < ptm_prot_num else ptm_prot_num
		outfile_name = "{}_temp_ptm_prots.txt".format(tax_id)
		outfile = "{}_{}_{}.txt".format(tax_id, str(i), str(endpos))
		outfile = open(outfile, 'w')
		cmd = 'curl -X GET "https://research.bioinformatics.udel.edu/iptmnet/api/browse?term_type=All&role=Substrate&organism={}&ptm_type=Phosphorylation&start_index={}&end_index={}" -H  "accept: application/json" -k > {}'.format(tax_id, str(i), str(endpos), outfile_name)
		system(cmd)
		prot_dicts = fileeval(outfile_name)
		remove(outfile_name)
		for prot_dict in prot_dicts:
			iptm_id = prot_dict["iptm_id"]
			uniprot_id = prot_dict["uniprot_ac"]
			cmd = 'curl -X GET "https://research.bioinformatics.udel.edu/iptmnet/api/{}/substrate" -H "accept: application/json" -k > temp.txt'.format(iptm_id)
			system(cmd)
			ptm_dicts = fileeval("temp.txt")
			outfile_write_strs(outfile, uniprot_id, iptm_id, ptm_dicts)
		outfile.close()

	print "UniProtID\tiPTMnetID\tposition\tresidue\tscore"
	for filename in listdir('./'):
		if filename.split('_')[0] == tax_id:
			for s in open(filename):
				print s.strip()


parser = OptionParser()
parser.add_option("-x", "--tax_id", help="TAX ID of species")
parser.add_option("-m", "--sample_size", help="Size of the single sample (default 100)", default = "100")
parser.add_option("-n", "--ptm_prot_num", help="Maximal number of proteins")
opt, args = parser.parse_args()

main(opt.tax_id, eval(opt.sample_size), eval(opt.ptm_prot_num))

