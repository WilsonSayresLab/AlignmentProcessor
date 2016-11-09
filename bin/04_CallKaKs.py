'''This program will run KaKs_Calculator on a directory. Must be in the 
AlignmentProcessor directory to run.


	Copyright 2016 by Shawn Rupp'''

import argparse
from subprocess import Popen
from shlex import split
from glob import glob
import os

DEVNULL = open(os.devnull, "w")

def calculateKaKs(indir, outdir, method):
	'''Calculates substition rates.'''
	print("\tRunning KaKs_Calculator...")
	files = glob(indir + "*.axt")
	for axt in files:
		with open(axt, "r") as infile:
			filename = axt.split("/")[-1]
			# Create output file
			outfile = (outdir + filename.split(".")[0] + ".kaks")
			ck = Popen(split("bin/KaKs_Calculator -i " + axt + " -o " +
							 outfile + " -m " + method), stdout = DEVNULL)
			ck.wait()
	return True

def compileKsKs(outdir):
	'''Prints Ka/Ks output as a single csv file.'''
	print("\tCompiling KaKs_Calculator output...")
	# Set counter so the header is only printed once
	count =  0
	files = glob(outdir + "*.kaks")
	# Create output file
	outfile = ""
	path = outdir.split("/")[:-2]
	for i in path:
		outfile += i + "/"
	outfile = outfile + "KaKs.csv"
	# Open input and output files
	with open(outfile, "w") as output:
		for kaks in files:
			with open(kaks, "r") as infile:
				filename = kaks.split("/")[-1]
				for line in infile:
					if count != 0:
						# Print data from remaining files
						if line.split("\t")[0] == "Sequence":
							pass
						else:
							output.write(filename.split(".")[0] + "\t" + 
										line.replace("\t",","))
					elif count == 0:
						#Print header from first file
						output.write("GeneID\t" + line.replace("\t",","))
						count += 1

def main():
	concatenate = False
	parser = argparse.ArgumentParser(description="This program will run \
KaKs_Calculator on a directory.")
	parser.add_argument("-i", help="Path to input file.")
	parser.add_argument("-o", help="Path to output file.")
	parser.add_argument("-m", help="Method for calculating Ka/Ks.") 
	# Parse arguments and assign to variables
	args = parser.parse_args()
	indir = args.i
	if indir[-1] != "/":
		indir += "/"
	outdir = args.o
	if outdir != "/":
		outdir += "/"
	method = args.m
	concatenate = calculateKaKs(indir, outdir, method)
	if concatenate == True:
		compileKsKs(outdir)

if __name__ == "__main__":
	main()
