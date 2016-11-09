'''This will take an aligned multiple FASTA file with multiple genes 
	and create individual FASTA alignment files for each gene and simplify the
	fasta headers.


	Copyright 2016 by Shawn Rupp'''

import argparse
from glob import glob

def splitFasta(infile, outdir):
	# Open input file and split into one alignment per gene
	print("\tSplitting fasta file into one file per gene...")
	passed = 0
	excluded = 0
	# Create log file
	log = ""
	path = outdir.split("/")[:-2]
	for i in path:
		log += i + "/"
	log += "splitFastaLog.txt"
	with open(log, "w") as logfile:
		logfile.write("Transcripts with only one species\n\n")
	# Parse input fasta
	with open(infile, "r") as fasta:
		newid = True
		seq = ""
		n = 0
		for line in fasta:
			if line != "\n":
				if line[0] == ">":
					line, geneid = convertHeader(line)
					# Concatenate lines for all species for each gene
					seq += str(line)
					# Determine number of sequences and species names
					n += 1
					if newid == True:
						# Set reference species ID as file name
						filename = geneid
						newid = False
				else:
					# Concatenate remaining lines
					seq += str(line)
			elif line == "\n" and newid == False:
				# Use empty lines to determine where genes end
				if n >= 2:
					# Print gene sequences to file if there are at least two
					# species and reset for next gene
					outfile = (outdir + filename + "." + str(n) + ".fa")
					with open(outfile, "w") as output:
						output.write(seq)
					newid = True
					seq = ""
					n = 0
					passed += 1
				elif n < 2:
					# Skip genes with only one sequence and save ID in log
					with open(log, "a") as logfile:
						logfile.write(geneid + "\n")
					excluded += 1
	with open(log, "a") as logfile:
		logfile.write(("\nTotal transcripts written to file: {}\n").format(passed))
		logfile.write(("Total transcripts with only one sequence: {}").format(excluded))

def convertHeader(line):
	'''Returns a header containing only the genome build name.'''
	if "_" in line:
		# Extract relevant data from UCSC header
		genebuild = line[1:].split()[0]
		genebuild = genebuild.split("_")
		build = ">" + str(genebuild[-1]) + "\n"
		geneid = str(genebuild[0].split(".")[0])
	else:
		# Extract build without ">" 
		build = ">" + line.split(".")[0][1:].rstrip() + "\n"
		geneid = str(line.split(".")[1]).rstrip()
	return build, geneid

def main():
	parser = argparse.ArgumentParser(description="This will take the \
aligned multiple FASTA file with multiple genes and create individual FASTA \
alignment files for each gene and simplify the fasta headers.")
	parser.add_argument("-i", help="path to input file.")
	parser.add_argument("-o", help="path to output directory.")
	args = parser.parse_args()
	infile = args.i
	outdir = args.o
	if outdir[-1] != "/":
		outdir += "/"
	splitFasta(infile, outdir)

if __name__ == "__main__":
	main()
