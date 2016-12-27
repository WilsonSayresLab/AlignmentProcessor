'''This program will run KaKs_Calculator on a directory. Must be in the 
AlignmentProcessor directory to run.


	Copyright 2016 by Shawn Rupp'''

from datetime import datetime
from sys import stdout
from glob import glob
from subprocess import Popen
from shlex import split
from functools import partial
from multiprocessing import Pool, cpu_count
import argparse
import os

# Define max number of threads and devnull for capturing stdout
MAXCPU = cpu_count()
DEVNULL = open(os.devnull, "w")

def calculateKaKs(indir, outdir, method, axt):
	'''Calculates substition rates.'''
	with open(axt, "r") as infile:
		filename = axt.split("/")[-1]
		# Create output file
		outfile = (outdir + filename.split(".")[0] + ".kaks")
		ck = Popen(split("bin/KaKs_Calculator -i " + axt + " -o " +
						 outfile + " -m " + method), stdout = DEVNULL)
		ck.wait()
	if ck.returncode() == 0:
		return True

def compileKsKs(outdir):
	'''Prints Ka/Ks output as a single tsv file.'''
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

#-----------------------------------------------------------------------------

def main():
	starttime = datetime.now()
	concatenate = False
	parser = argparse.ArgumentParser(description="This program will run \
KaKs_Calculator on a directory.")
	parser.add_argument("-i", help = "Path to input file.")
	parser.add_argument("-o", help = "Path to output file.")
	parser.add_argument("-m", help = "Method for calculating Ka/Ks.") 
	parser.add_argument("-t", type = int, help = "Number of threads.")
	# Parse arguments and assign to variables
	args = parser.parse_args()
	indir = args.i
	if indir[-1] != "/":
		indir += "/"
	outdir = args.o
	if outdir != "/":
		outdir += "/"
	method = args.m
	cpu = args.t
	if cpu > MAXCPU:
		cpu = MAXCPU
	# Call Ka/Ks_Calculator in parallel.
	genes = glob(indir + "*.axt")
	l = int(len(genes))
	pool = Pool(processes = cpu)
	func = partial(calculateKaKs, indir, outdir, method)
	print("\tRunning KaKs_Caclulator with", str(cpu), "threads....")
	rcml = pool.imap(func, genes)
	pool.close()
	pool.join()
	# Compile output
	compileKsKs(outdir)
	print("\tKaKs_Calculator runtime: ", datetime.now() - starttime)

if __name__ == "__main__":
	main()
