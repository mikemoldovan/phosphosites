"""
filter_by_score

Filter phosphosite TSV file by score

Input:
1. TSV phosphosite file
2. Score

Output:
TSV filtered by score 
"""

from optparse import OptionParser

def main(tsv_filename, min_score):
	with open(tsv_filename) as tsv:
		for s in tsv:
			s = s.strip().split()
			if s[0] == "UniProtID":
				print '\t'.join(s)
				continue
			if eval(s[4]) >= min_score:
				print '\t'.join(s)


parser = OptionParser()
parser.add_option("-i", "--tsv_filename", help="TSV phosphosite file")
parser.add_option("-s", "--min_score", help="Score")
opt, args = parser.parse_args()


main(opt.tsv_filename, eval(opt.min_score))