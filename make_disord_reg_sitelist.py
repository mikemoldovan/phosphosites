"""
make_disord_reg_sitelist

Infer sites from ordered/disordered regions

Input:
1. File containing a list of phosphosites
2. File containing disordered regions
3. Flag: d for disordered and o for ordered sites

Output:
substitution matrix
"""


from optparse import OptionParser
from disord_reg_adlib import *
from sys import stderr

def main(psite_filename, dr_filename):
#	print "UniProtID\tiPTMnetID\tposition\tresidue\tscore"
	psite_dict = make_psite_dict(psite_filename)
	dr_dict = build_dr_dict(dr_filename)
	with open(psite_filename) as psites:
		for s in psites:
			s = s.strip().split()
			if s[0] == "UniProtID":
				continue
			try:
				disord_reg_crds = dr_dict[s[0]]
			except KeyError:
				pass
#				stderr.write(s[0]+'\n')
			crd = eval(s[2]) - 1
			for crd_pair in disord_reg_crds:
				if crd >= crd_pair[0] and crd < crd_pair[1]:
					print '\t'.join(s)



parser = OptionParser()
parser.add_option("-p", "--psite_filename", help="File containing a list of phosphosites")
parser.add_option("-d", "--disord_regs_filename", help="File containing disordered (ordered) regions")
opt, args = parser.parse_args()


main(opt.psite_filename, opt.disord_regs_filename)