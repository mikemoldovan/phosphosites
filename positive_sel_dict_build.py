"""
positive_sel_dict_build

Builds a dictionary with numbers of sites unser positive selection

Input:
1. Phosphosite file
2. Disordered region file
3. Directory with alignments
4. Directory with GODON outputs
5. Reference species ID
"""

from os import listdir
from optparse import OptionParser
from Bio import SeqIO
from disord_reg_adlib import build_dr_dict, make_psite_dict2


genecode = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
            "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
            "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
            "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
            "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
            "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
            "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
            "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
            "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
            "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
            "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
            "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
            "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
            "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
            "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
            "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}


def make_possel_dict(godon_outfile):
	possel_dict = dict()
	with open(godon_outfile) as goout:
		for s in goout:
			if s[0].isdigit():
				s = s.strip().split()
				possel_dict[eval(s[0]) - 1] = s
	return possel_dict


def make_aldict(alfile, refspec_id):
	aldict = dict()
	refspec_seq = ""
	with open(alfile) as inhandle:
		for record in SeqIO.parse(inhandle, "fasta"):
			s = ""
			for c in record.seq:
				s += c
			aldict[record.id] = s
			if record.id == refspec_id:
				refspec_seq = s
	return aldict, refspec_seq


def get_col_aacids(al_crd, aldict):
	global genecode
	col_aacids = []
	for k in aldict.keys():
		codon = aldict[k][al_crd : al_crd + 3]
		if 'X' in codon:
			continue
		if '-' in codon:
			continue
		aacid = genecode[codon]
		if aacid == 's':
			print aacid
		if aacid not in col_aacids:
			col_aacids.append(aacid)
	return ''.join(sorted(col_aacids))

def possel_append(possel_dict,
			      aldict,
			      refspec_seq,
			      prot_id,
			      refspec_id,
			      possel_numdict,
			      dr_dict,
			      psite_dict):
	global genecode
	count_al = -1
	count_seq = -1
	current_dr = 0
	for i in range(0, len(refspec_seq), 3):
		count_al += 1
		if '-' not in refspec_seq[i:i+3]:
			count_seq += 1
		else:
			continue
		if 'X' in refspec_seq[i:i+3]:
			continue
		if not dr_dict[prot_id]:
			break
#		try:
#			a = count_seq > dr_dict[prot_id][current_dr][1]
#		except:
#			print(count_seq, dr_dict[prot_id])
		if count_seq > dr_dict[prot_id][current_dr][1]:
			current_dr += 1
			if current_dr == len(dr_dict[prot_id]):
				break
		if count_seq >= dr_dict[prot_id][current_dr][0] and count_seq < dr_dict[prot_id][current_dr][1]:
			if not psite_dict.get(prot_id):
				continue
			if psite_dict[prot_id].get(count_seq):
				ref_aacid = 'p' + genecode[refspec_seq[i:i+3]]
			else:
				ref_aacid = genecode[refspec_seq[i:i+3]]
			if not possel_dict.get(count_al):
				continue
			
			col_aacids = get_col_aacids(i, aldict)
			if possel_numdict.get(ref_aacid):
				if possel_numdict[ref_aacid].get(col_aacids):
					possel_numdict[ref_aacid][col_aacids] += 1
				else:
					possel_numdict[ref_aacid][col_aacids] = 1
			else:
				possel_numdict[ref_aacid] = {col_aacids : 1}


def main(psite_filename, dr_filename, aldir_name, godon_out_dir, refspec_id):
	psite_dict = make_psite_dict2(psite_filename)
	dr_dict = build_dr_dict(dr_filename)
	possel_numdict = dict()
	for godon_outfile in listdir(godon_out_dir):
		possel_dict = make_possel_dict(godon_out_dir + godon_outfile)
		if len(possel_dict.keys()) == 0:
			continue
		alfile_name = aldir_name + godon_outfile.split('.')[0] + ".fa"
		aldict, refspec_seq = make_aldict(alfile_name, refspec_id)
		if not refspec_seq:
			continue
		prot_id = godon_outfile.split('.')[0][:-4]
		possel_append(possel_dict,
			          aldict,
			          refspec_seq,
			          prot_id,
			          refspec_id,
			          possel_numdict,
			          dr_dict,
			          psite_dict)
	print(possel_numdict)
	print(list(possel_numdict.keys()))


parser = OptionParser()
parser.add_option("-p", "--psite_filename", help="Phosphosite file")
parser.add_option("-d", "--dr_filename", help="Disordered region file")
parser.add_option("-a", "--aldir_name", help="Directory with alignments")
parser.add_option("-g", "--godon_out_dir", help="Directory with GODON outputs")
parser.add_option("-r", "--refspec_id", help="Reference species ID")
opt, args = parser.parse_args()


main(opt.psite_filename, opt.dr_filename, opt.aldir_name, opt.godon_out_dir, opt.refspec_id)