"""
obtain_ord_reg_list

Make ordered region file from the file with disordered regions

Input:
1. File with disordered regions
2. Proteome FASTA

Output:
ordered regions
"""
from optparse import OptionParser
from disord_reg_adlib import *

def main(dr_filename, proteome_fasta):
	fastadict = make_fastadict(proteome_fasta)
	with open(dr_filename) as dr:
		for protrec in dr:
			protrec = protrec.strip().split()
			protlen = len(fastadict[protrec[0]])
			if len(protrec) < 2:
				print "{}\t{}-{}".format(protrec[0], 1, protlen)
				continue
			crd_str = ""
			crd_list = protrec[1].split('|')
			crd_0 = crd_list[0].split('-')[0]
			if eval(crd_0) > 2:
				crd_str += "{}-{}|".format(1, eval(crd_0) - 1)
			crds1 = crd_list[0].split('-')
			for i in range(1, len(crd_list)):
				crds0 = crd_list[i-1].split('-')
				crds1 = crd_list[i].split('-')
				crd_str += "{}-{}|".format(eval(crds0[1]) + 1, eval(crds1[0]) - 1)
			if eval(crds1[1]) < (protlen - 1):
				crd_str += "{}-{}|".format(eval(crds1[1]) + 1, protlen)
			print "{}\t{}".format(protrec[0], crd_str[:-1])


parser = OptionParser()
parser.add_option("-d", "--proteome_fasta", help="Proteome file")
parser.add_option("-r", "--dr_filename", help="Disordered region file")
opt, args = parser.parse_args()


main(opt.dr_filename, opt.proteome_fasta)