#!/usr/bin/env python3

# NOTE: this file was modified to fit the SQANTI-SIM pipeline
# Contributor: Jorge Mestre

import sys,time,argparse


def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	generate_expr_matrix(args.input,args.input_expr,args.output, args.number)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def generate_expr_matrix(input_gpd_fl,input_txt_fl,output_expr_mtx, n_reads):
	mill_reads = n_reads/1000000
	# parse txt
	dic_iso_expr = {}
	skip = input_txt_fl.readline() # header
	for line in input_txt_fl:
		values = line.strip().split() 
		iso_id = values[0]
		expr_v = values[2] # expr_v is TPM now
		dic_iso_expr[iso_id] = int(round(mill_reads * float(expr_v)))
	
	for line in input_gpd_fl:
		iso_id = line.strip().split()[1]
		if iso_id in dic_iso_expr:
			output_expr_mtx.write(line.strip() + "\t" + str(dic_iso_expr[iso_id]) + '\n')

	input_txt_fl.close()
	input_gpd_fl.close()
	output_expr_mtx.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Randomly generate read count for each isoform based on negative binomial (NB) distribution. Read count is shown in last column of output file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-e','--input_expr',type=argparse.FileType('r'),required=True,help="Input: expression txt file (first colunm is isoform ID, second colunmn is expression value)")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd + read count file")
	parser.add_argument('--number', '-n', help='number of reads to generate', default=20000, type=int)
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
