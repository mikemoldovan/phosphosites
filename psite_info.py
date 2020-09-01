"""
psite_info
get information about a set of phosphosites

input:
file in format: niprotID\IPTMnetID\position\residue\score (iptmnet_info output)

output:
1. Number of phosphorilated proteins
2. Total number of phosphosites
3. Numbers of phosphosites for each score
4. Numbers of phosphosites for each aminoacid for each score
"""


from optparse import OptionParser

def main(psite_info_filename):
	score_dict = dict()
	by_aa_score_dict = dict()
	phosphoprot_dict = dict()
	phosphoprot_dict_2 = dict()
	total_psites = 0
	with open(psite_info_filename) as psite_file:
		for s in psite_file:
			s = s.strip().split()
			if s[0] == "UniProtID":
				continue

			total_psites += 1

			if score_dict.get(s[4]):
				score_dict[s[4]] += 1
			else:
				score_dict[s[4]] = 1

			phosphoprot_dict[s[0]] = True

			if eval(s[4]) > 1:
				phosphoprot_dict_2[s[0]] = True

			if by_aa_score_dict.get(s[3]):
				if by_aa_score_dict[s[3]].get(s[4]):
					by_aa_score_dict[s[3]][s[4]] += 1
				else:
					by_aa_score_dict[s[3]][s[4]] = 1
			else:
				by_aa_score_dict[s[3]] = dict()
				by_aa_score_dict[s[3]][s[4]] = 1

	print "Total psites:\t{}".format(total_psites)
	print "Phosphorolated proteins:\t{}".format(len(phosphoprot_dict.keys()))
	print "Good phosphorolated proteins:\t{}".format(len(phosphoprot_dict_2.keys()))
	for score in sorted(score_dict.keys()):
		print "Score {}:\t{}".format(score, score_dict[score])
	for aa in by_aa_score_dict.keys():
		for score in sorted(by_aa_score_dict[aa].keys()):
			print "{} with score {}\t{}".format(aa, score, by_aa_score_dict[aa][score])


parser = OptionParser()
parser.add_option("-i", "--psite_info_filename", help="file in format: niprotID\IPTMnetID\position\residue\score (iptmnet_info output)")
opt, args = parser.parse_args()

main(opt.psite_info_filename)





