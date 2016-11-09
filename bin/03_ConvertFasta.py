'''This program will convert single fasta alignemnt files to either axt or 
sequential phylip format.


	Copyright 2016 by Shawn Rupp'''

import argparse
from glob import glob
import os

def convert(indir, outdir, axt, phylip):
	'''Open all input files in the directory and call conversion function'''
	files = glob(indir + "*")
	if axt == True:
		fastaToAxt(files, outdir)
	elif phylip == True:
		fastaToPhylip(files, outdir)

def fastaToAxt(files, outdir):
	'''Reformats fasta alignment to axt format.'''
	print("\tConverting fasta files to axt...")
	for fasta in files:
		# Create output file:
		filename = fasta.split("/")[-1]
		geneid = filename.split(".")[0]
		n = filename.split(".")[1]
		outfile = outdir + geneid + ".axt"
		with open(fasta, "r") as infile:
			axt = []
			header1 = ""
			seq1 = ""
			# Parse input file
			for line in infile:
				if line[0] == ">":
					if not header1:
						header1 = line[1:].strip()
					elif header1:
						header2 = line[1:].strip()
				elif line[0] != ">":
					if not seq1:
						seq1 = line.strip()
					elif seq1:
						seq2 = line.strip()
		# Assemble axt header and add sequences
		header = header1 + "-" + header2
		axt.append(header)
		axt.append(seq1)
		axt.append(seq2)
		with open(outfile, "w") as output:
				for line in axt:
					output.write(line + "\n")

def fastaToPhylip(files, outdir):
	'''Reformats fasta alignment to sequential phylip format.'''
	print("\tConverting fasta files to phylip...")
	for fasta in files:
		seqs = []
		n = 0
		filename = fasta.split("/")[-1]
		geneid = filename.split(".")[0]
		with open(fasta, "r") as infile:
			# Parse fasta file and convert to phylip
			for line in infile:
				if line[0] == ">":
					species = line[1:].strip()
					if len(species) > 10:
						species = species[:11]
					n += 1
				else:
					if n == 1:
						l = len(line.strip())
					seqs.append(species + "  " + line.strip())
		n = str(n)
		outfile = outdir + geneid + "." + str(n) + ".phylip"
		with open(outfile, "w") as output:
			output.write(" " + n + " " + str(l) + "\n")
			for i in seqs:
				output.write(i + "\n")

def main():
	parser = argparse.ArgumentParser(description="This program will convert \
single fasta alignemnt files to either axt or sequential phylip format.")
	parser.add_argument("-i", help="Path to input file.")
	parser.add_argument("-o", help="Path to output directory.")
	parser.add_argument("--axt", action="store_true", 
help="Convert files to axt format.")
	parser.add_argument("--phylip", action="store_true",
help="Convert files to sequential phylip format.")
	# Parse arguments and assign to variables
	args = parser.parse_args()
	indir = args.i
	if indir[-1] != "/":
		indir += "/"
	outdir = args.o
	if outdir != "/":
		outdir += "/"
	axt = args.axt
	phylip = args.phylip
	if axt == True and phylip == True:
		print("\tPlease specify only one file format to convert to.")
		quit()
	convert(indir, outdir, axt, phylip)

if __name__ == "__main__":
	main()
