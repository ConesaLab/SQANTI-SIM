#!/usr/bin/env python3

# NOTE: this file was modified to fit the SQANTI-SIM pipeline
# Contributor: Jorge Mestre

import sys,time,argparse
import numpy as np


def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	convert_gpd_to_fasta(args.input_fasta,args.input_gpd,args.output_fasta)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def convert_gpd_to_fasta(input_fa_fl,input_gpd_fl,output_fa_fl):
	# --- parse genome fasta file ---
	dic_chr_seq = {}
	chr = ""
	seq_list = [""]
	for line in input_fa_fl:
		if line.startswith(">"):
			dic_chr_seq[chr] = "".join(seq_list)
			seq_list = []
			chr = line.strip().split(" ")[0][1:]
		else:
			seq_list.append(line.strip())
	dic_chr_seq[chr] = "".join(seq_list) # add last chromosome
	del dic_chr_seq[""] # delete null key/value pair

	# --- parse gpd annotation and generate fasta for each isoform ---
	dic_gatc_com = {"G":"C","A":"T","T":"A","C":"G"}
	dic_nogatc_code = {"U":["T"],"I":["G"],"R":["A","G"],"Y":["C","T"],"S":["G","C"],"W":["A","T"],"K":["G","T"],"M":["A","C"],"B":["C","G","T"],"D":["A","G","T"],"H":["A","C","T"],"V":["A","C","G"],"N":["A","C","G","T"],"X":["A","C","G","T"]}
	for line in input_gpd_fl:
		gene,iso,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.strip().split("\t")[:11]
		seq_list = []
		for i in range(0,int(exon_number)):
			start = int(exon_start.split(",")[i])
			end = int(exon_end.split(",")[i])
			seq_list.append(dic_chr_seq[chr][start:end])
		seq = "".join(seq_list)
		seq = seq.upper()

		up_seq_list = list(seq) # non AGCT nucleotide
		for i in range(len(up_seq_list)):
			if up_seq_list[i] in dic_nogatc_code:
				up_seq_list[i] = np.random.choice(dic_nogatc_code[up_seq_list[i]])

		if strand == "-": # minus strand, get reverse complementary sequence
			seq_com = ""
			for nuc in up_seq_list:
				seq_com += dic_gatc_com[nuc]
			seq = seq_com[::-1]
		else:
			seq = "".join(up_seq_list)
		output_fa_fl.write(">" + iso + '\n')
		output_fa_fl.write(seq + '\n')

	input_fa_fl.close()
	input_gpd_fl.close()
	output_fa_fl.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Generate fasta sequence for each isoform.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-a','--input_fasta',type=argparse.FileType('r'),required=True,help="Input: genome fasta file")
	parser.add_argument('-g','--input_gpd',type=argparse.FileType('r'),required=True,help="Input: annotation gpd file")
	parser.add_argument('-o','--output_fasta',type=argparse.FileType('w'),required=True,help="Output: isoform fasta file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
