"""
build_and_trim_alignment

Takes a protein sequence file and a corresponding nucleotide sequence file and
performs several stages of alignment building and refinement:

1. Discard all nucleotide sequences which do not correspond to protein sequences
2. Trim stop codons
3. Substitute all codons with 'X' with 'XXX'
4. Discard files with <6 sequences
5. Build alignment
6. Discard all columns with <4 residues
7. Compute Guidance scores for the alignment
8. Substitute all residues with <0.93 score by 'X'
"""

from Bio import SeqIO
from Bio.Seq import Seq
from optparse import OptionParser
from os import system

######## Discard all nucleotide sequences which do not correspond to protein sequences #######

def make_seqdict(fafile):
	seqdict = dict()
	with open(fafile) as infasta:
		for record in SeqIO.parse(infasta, "fasta"):
			seqdict[record.id] = record
	return seqdict

def corresp(protseq, nuclseq):
	if len(nuclseq) < 30:
		return False
	if len(protseq) < 10:
		return False
	if 'X' not in nuclseq[:30]:
		s_n = ""
		for c in nuclseq[:30].translate():
			s_n += c
		s_p = ""
		for c in protseq[:10]:
			s_p += c
		if s_n == s_p:
			return True
		else:
			return False
	else:
		corresp(protseq[10:], nuclseq[30:])

def discard_nucl_seqs(protseq_fasta, nucl_seq_fasta, outfile_name):
	prot_seqdict = make_seqdict(protseq_fasta)
	nucl_seqdict = make_seqdict(nucl_seq_fasta)
	with open(outfile_name, 'w') as outfile:
		for seq_id in nucl_seqdict.keys():
			protseq = prot_seqdict[seq_id].seq
			nuclseq = nucl_seqdict[seq_id].seq
			if corresp(protseq, nuclseq):
#				print(corresp(protseq, nuclseq))
				SeqIO.write(nucl_seqdict[seq_id], outfile, "fasta")

##### Trim stop codons #####

def discard_stop_codon(seq, stop_codon):
	for i in range(0, len(seq), 3):
		part_seq = seq[i] + seq[i + 1] + seq[i + 2]
		if part_seq == stop_codon:
			return seq[:i]
	return seq
	

def trim_stop_codons(nucl_seq_fasta, outfile_name):
	with open(nucl_seq_fasta) as infasta, open(outfile_name, 'w') as outfile:
		for record in SeqIO.parse(nucl_seq_fasta, "fasta"):
			newseq = record.seq
			while 'X' in newseq[-3:]:
				newseq = newseq[:-3]
			newseq = discard_stop_codon(newseq, "TGA")
			newseq = discard_stop_codon(newseq, "TAA")
			newseq = discard_stop_codon(newseq, "TAG")
			record.seq = newseq
			SeqIO.write(record, outfile, "fasta")


##### Substitute all codons with 'X' with 'XXX' #####
##### Discard files with <6 sequences #####

def align_discard(nucl_seq_fasta, thresh_number):
	count = 0
	with open(nucl_seq_fasta) as infasta:
		for s in infasta:
			if s.startswith('>'):
				count += 1
	if count < thresh_number:
		return True
	return False


##### Build alignment #####

def build_align(nucl_seq_fasta, outfile_name):
	system("perl ~/tools/translatorx.pl -i tmp2.fa -o tmp.fa")


##### Discard all columns with <4 residues #####

def make_seqdict2(fafile):
	seqdict = dict()
	with open(fafile) as infasta:
		for record in SeqIO.parse(infasta, "fasta"):
			s = ""
			for c in record.seq:
				s += c
			seqdict[record.id] = list(s)
	return seqdict

def small_column_discard(nucl_al_fasta, outfile_name):
	seqdict = make_seqdict2(nucl_al_fasta)
	keylist = list(seqdict.keys())
	align_len = len(seqdict[keylist[0]])
	for i in range(align_len):
		resnum = 0
		for key in keylist:
			if seqdict[key][i] != '-':
				resnum += 1
		if resnum < 4:
			for key in keylist:
				if seqdict[key][i] != '-':
					seqdict[key][i] = 'X'
	with open(outfile_name, 'w') as outhandle:
		for key in keylist:
			outhandle.write('>' + key + '\n')
			outhandle.write(''.join(seqdict[key]) + '\n')


def main(protseq_fasta, nucl_seq_fasta, outfile_name):
	discard_nucl_seqs(protseq_fasta, nucl_seq_fasta, "tmp1.fa")
	trim_stop_codons("tmp1.fa", "tmp2.fa")
	if align_discard("tmp2.fa", 6):
		system("rm tmp*")
		return 0
	build_align("tmp2.fa", outfile_name)
	system("mv tmp.fa.nt_ali.fasta tmp3.fa")
	small_column_discard("tmp3.fa", outfile_name)
	system("rm tmp*")


parser = OptionParser()
parser.add_option("-n", "--nucl_seq_fasta", help="Nucleotide sequence file")
parser.add_option("-p", "--protseq_fasta", help="Corresponding protein sequence file")
parser.add_option("-o", "--outfile_name", help="Name of the output file")
opt, args = parser.parse_args()

main(opt.protseq_fasta, opt.nucl_seq_fasta, opt.outfile_name)
