#from /home/mmoldovan/phosphosites/scripts/disord_reg_adlib import *
from Bio import SeqIO

fastadict1 = dict()
fastadict2 = dict()
iddict = dict()

with open("proteomes/mouse_proteome.fa") as inhandle:
	for record in SeqIO.parse(inhandle, "fasta"):
		s = ""
		for c in record.seq:
			s += c
		fastadict1[record.id.split('|')[2]] = s
		fastadict2[record.id.split('|')[1]] = s
		iddict[record.id.split('|')[1]] = record.id.split('|')[2]


with open("mouse_tissues_raw.txt") as inhandle:
	for psite_record in inhandle:
		_str = psite_record
	for psite_record in _str.split("IPI:IPI"):
		s = psite_record[:-1]
		psite_record = s.strip().split('\t')
		infadict = False
		if "SWISS-PROT" not in psite_record[0]:
			continue
		for k in fastadict2.keys():
			if k in psite_record[0]:
				print iddict[k] + '\t' + k + '\t' + psite_record[2] + '\t' + psite_record[1] + '\t' + '\t'.join(psite_record[3:])
				infadict = True
				break
#		if not infadict:
#			print "*" + s


