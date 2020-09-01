"""
pondr_run_2

runs PONDR VSL and yields a file with disordered regions

Input:
1. PAML output directory
2. Job ID

Output:
File in format: uniprot_id\tseq_id_in_paml_run\tn1-n2;n3-n4;...
n1-n2 -- disordered region
"""


from os import system, remove, listdir
from optparse import OptionParser
from Bio import SeqIO


def print_disordered(seq, uniprot_id, seq_id_in_paml_run, job_id):
	with open("temp_{}.txt".format(job_id), 'w') as temp:
		temp.write(seq)
	t = open("temp_out_{}.txt".format(job_id), 'w')
	t.close()
	system("java -jar ~/tools/VSL2/VSL2.jar -s:temp_{}.txt > temp_out_{}.txt".format(job_id, job_id))
	flag = False
	outstr = uniprot_id + '\t' + seq_id_in_paml_run + '\t'
	disordered_regs = []
	with open("temp_out_{}.txt".format(job_id)) as outfile:
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
	remove("temp_{}.txt".format(job_id))
	remove("temp_out_{}.txt".format(job_id))


def main(paml_outdir, job_id):
	for filename in listdir(paml_outdir):
		if not filename.endswith(".fa"):
			continue
		uniprot_id = filename[:-3]
		proteome_fasta = paml_outdir + filename
		with open(proteome_fasta) as infasta:
			for record in SeqIO.parse(infasta, "fasta"):
				seq_id_in_paml_run = record.id
				seq = ""
				for c in record.seq:
					if c != '-':
						if c == 'X':
							seq += 'A'
						else:
							seq += c
				print_disordered(seq, uniprot_id, seq_id_in_paml_run, job_id)


parser = OptionParser()
parser.add_option("-i", "--paml_outdir", help="PAML output directory")
parser.add_option("-j", "--job_id", help="Job ID", default='0')
opt, args = parser.parse_args()


main(opt.paml_outdir, opt.job_id)
