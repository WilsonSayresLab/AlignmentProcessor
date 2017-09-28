'''This program will mask stop codons.

	Copyright 2016 by Shawn Rupp'''

import argparse
from glob import glob
from collections import OrderedDict
import os

def openFiles(indir, outdir):
	'''Opens all input files in the directory, performs filtering steps, 
and writes to file.'''
	print("\tMasking stop codons...")
	# Iterate through files
	files = glob(indir + "*")
	for fasta in files:
		filename = fasta.split("/")[-1]
		geneid = filename.split(".")[0]
		# Extract number of sequnces from file name
		n = int(filename.split(".")[1])
		seqs = seqDict(fasta)
		seqs = removeStops(seqs)
		# Create output file
		outfile = (outdir + geneid + "." + str(n) + ".filtered.fa")
		writeDict(outfile, seqs)

def seqDict(fasta):
	'''Convert fasta into separate sequence objects, determine sequence names
	and create dictionary entries for each set of codons'''
	seqs = OrderedDict()
	with open(fasta, "r") as infile:
		for line in infile:
			if line[0] == ">":
				species = line[1:].strip()
			if line[0] != ">":
				codons = []
				seq = line.strip()
				for i in range(0, len(seq), 3):
					codons.append(seq[i:i +3])
					i += 3
				seqs[species] = codons
	return seqs

def removeStops(seqs):
	'''This will replace internal stop codons with gaps.'''
	for species in seqs:
		for idx,codon in enumerate(seqs[species]):
			if codon == "TAA" or codon == "TAG" or codon == "TGA":
					# Replace stops with dashes
					seqs[species][idx] = "---"
	return seqs
	 
def writeDict(outfile, seqs):
	# Convert dictionary values to string and write to file
	with open(outfile, "w") as output:	 
		for species in seqs:
			output.write(">" + str(species) + "\n")
			seq = ""
			for codon in seqs[species]:
				seq += str(codon)
			output.write(seq + "\n")

def main():
	parser = argparse.ArgumentParser(description="This script will mask stop \
codons in a fasta alignment.")
	parser.add_argument("-i", help="Path to input directory.")
	parser.add_argument("-o", help="Path to output directory.")
	# Parse arguments and assign to variables
	args = parser.parse_args()
	indir = args.i
	if indir[-1] != "/":
		indir += "/"
	outdir = args.o
	if outdir[-1] != "/":
		outdir += "/"
	openFiles(indir, outdir)

if __name__ == "__main__":
	main()
