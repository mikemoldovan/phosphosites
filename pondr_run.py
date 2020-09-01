"""
pondr_run

runs PONDR VSL and yields a file with disordered regions

Input:
FASTA file with proteome

Output:
File in format: seq_id\tn1-n2;n3-n4;...
n1-n2 -- disordered region
"""


from os import system, remove
from optparse import OptionParser
from Bio import SeqIO


def print_disordered(seq, uniprot_id):
	with open("temp.txt", 'w') as temp:
		temp.write(seq)
	system("java -jar ~/tools/VSL2/VSL2.jar -s:temp.txt > temp_out.txt")
	flag = False
	outstr = uniprot_id + '\t'
	disordered_regs = []
	with open("temp_out.txt") as outfile:
		for s in outfile:
			if "Predicted Disordered Regions" in s:
				flag = True
				continue
			if flag and (not s[0].isdigit()):
				flag = False
				break
			if flag:
				s = s.strip()
				disordered_regs.append(s)
	print outstr + '|'.join(disordered_regs)
	remove("temp.txt")
	remove("temp_out.txt")


def main(proteome_fasta):
	with open(proteome_fasta) as infasta:
		for record in SeqIO.parse(infasta, "fasta"):
			uniprot_id = record.id.split('|')[2]
			seq = ""
			for c in record.seq:
				seq += c
			print_disordered(seq, uniprot_id)


parser = OptionParser()
parser.add_option("-i", "--proteome_fasta", help="FASTA file with proteome")
opt, args = parser.parse_args()


main(opt.proteome_fasta)
