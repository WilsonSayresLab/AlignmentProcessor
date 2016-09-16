'''This program executes parseFastaIntoAXT.pl on an entire direcotry,
allowing all of the contents of the directory to be converted to axt files.


	Copyright 2016 by Shawn Rupp'''

from sys import argv
from subprocess import Popen
from shlex import split
from glob import glob
import os

def fastToAxt(path):
	'''Open all input files in the directory and convert to axt script for
 use with KaKs_Calculator'''
	inpath = path + "05_ReplaceStopCodons/" + "*.rmStops"
	files = glob(inpath)
	for file in files:
		axt = convert(file)
		# Create output file:
		filename = file.split("/")[-1]
		geneid = filename.split(".")[0]
		outfile = (path + "06_axtFiles/" + geneid + ".axt")
		with open(outfile, "w") as output:
			for line in axt:
				output.write(line + "\n")

def convert(file):
	'''Reformats fasta alignment to axt format.'''
	with open(file, "r") as infile:
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
	return axt

def main():
	if argv[1] == "-h" or argv[1] == "--help":
		print("Usage: python 06_FASTAtoAXT.py \
<path to inut and output directories>")
		quit()
	else:
		path = argv[1]
		fastToAxt(path)

if __name__ == "__main__":
	main()
