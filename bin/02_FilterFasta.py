'''This program will ensure that the open reading frame of a fasta alignment
	is conserved in each species relative to a reference species. It will 
	ensure that each sequence contains a specied threshold of nucleotides 
	before removing stop codons.

	Copyright 2016 by Shawn Rupp'''

import argparse
from glob import glob
from collections import OrderedDict
import os

def openFiles(indir, outdir, ref, percent, retainStops):
	'''Opens all input files in the directory, performs filtering steps, 
and writes to file.'''
	print("\tFiltering fasta alignments...")
	printed = 0
	removed = 0
	outfiles = glob(outdir + "*")
	if len(outfiles) > 0:
		# Delete contents of folder to avoid index errors
		for i in outfiles:
			os.remove(i)
	# Create log file
	log = ""
	path = outdir.split("/")[:-2]
	for i in path:
		log += i + "/"
	log += "filterLog.txt"
	with open(log, "w") as logfile:
		logfile.write("Transcript ID\tSpecies\tReason for Removal\n")
	# Iterate through files
	files = glob(indir + "*")
	for fasta in files:
		check = False
		fix = False
		count = False
		stops = False
		save = False
		filename = fasta.split("/")[-1]
		geneid = filename.split(".")[0]
		# Extract number of sequnces from file name
		n = int(filename.split(".")[1])
		check, seqs = seqDict(fasta, n)
		if check == True:
			try:
				fix, newseq = checkFrame(seqs, ref)
			except KeyError:
				# Skip genes which don't contain the reference species
				pass
		if fix == True:
			stops, newseq = fixCodons(newseq)
		if stops == True:
			try:
				count, newseq, n = removeStops(n, newseq, retainStops, 
												geneid, log)
			except TypeError:
				# Alignments with one species raise NoneType not iterable
				with open(log, "a") as logfile:
					logfile.write(geneid + "\tAll\tOnly contains 1 species \n")
		if count == True:
			try:
				save, newseq, n = countBases(n, newseq, percent, geneid, log)
			except TypeError:
				# Alignments with one species raise NoneType not iterable
				with open(log, "a") as logfile:
					logfile.write(geneid + "\tAll\tOnly contains 1 species \n")
		if save == True:
			# Create output file with updated n:
			outfile = (outdir + geneid + "." + str(n) + ".filtered.fa")
			writeDict(outfile, newseq)
		else:
			pass

def seqDict(fasta, n):
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
		return True, seqs

#-----------------------------------------------------------------------------

def checkFrame(seqs, ref):
	'''Removes gaps in the reference sequence introduced by
	alignment and removes corresponding sites in other species'''
	for codon in seqs[ref]: 
		if codon == "---":
			pass
		else:
			if "-" in str(codon):
				# Remove gaps from reference sequence
				i = codon.index("-")
				idx = seqs[ref].index(codon)
				codon = codon[:i] + codon[i + 1:]
				for key in seqs:
					# Remove corresponding gaps in other sequences
					triplet = seqs[key][idx]
					triplet = triplet[:i] + triplet[i + 1:]
	return True, seqs
						
def fixCodons(newseq):
	'''Replaces codons with missing nucleotides with gaps to 
remove unknown amino acids'''
	for species in  newseq:
		for idx, codon in enumerate(newseq[species]):
			if codon == "---":
				# Skip unknown codons
				pass
			else:
				for i,char in enumerate(codon):
					if char == "A" or char == "T" or char == "C" or char == "G":
						pass
					elif char == "-":
						# Replace ambiguous codons
						newseq[species][idx] = "---"
						break
	return True, newseq

def removeStops(n, newseq, retainStops, geneid, log):
	'''This will replace internal stop codons with gaps, create a list
of genes with internal stop codons, and remove those sequences.'''
	remove = []
	for species in newseq:
		# Replace terminal stop codons so the program can identify
		# remaining internal stops
		record = True
		try:
			if (newseq[species][-1] == "TAA" or newseq[species][-1] == "TAG" or
				newseq[species][-1] == "TGA"):
				newseq[species][-1] = "---"
		except (IndexError, ValueError) as e:
			# Remove empty sequences
			n -= 1
			del newseq[species] 
			with open(log, "a") as logfile:
				logfile.write(geneid + "\t" + species + 
							"\tLow Nucleotide Count\n")
			break
		for idx,codon in enumerate(newseq[species]):
			if codon == "TAA" or codon == "TAG" or codon == "TGA":
				# Record sequence once if it contains any stops 
				if record == True:
					record = False
					remove.append(species)
					with open(log, "a") as logfile:
						logfile.write(geneid + "\t" + species + 
									"\tPremature Stop Codon\n")
				if retainStops == True:
					# Replace stops with dashes
					newseq[species][idx] = "---"
	if retainStops == False:
		# Remove key after iterating to avoid using multiple breaks
		for species in remove:
			n -= 1
			del newseq[species]
	if n >= 2:
		return True, newseq, n
	elif n < 2:
		# Skip files with only one sequence left
		pass

def countBases(n, newseq, percent, geneid, log):
	'''Counts the number of nucleotides and only writes the sequence to an
output file if they compose greater than the cutoff threshold of the sequence'''
	failed = []
	for species in newseq:
		seq = ""
		for i in newseq[species]:
			seq += str(i)
		# Get count of all nucleotides
		gcount = seq.count("G")
		ccount = seq.count("C")
		acount = seq.count("A")
		tcount = seq.count("T")
		aligned = acount + tcount + ccount + gcount
		try:
			# Determine whether or not the sequence passes the treshold
			if aligned/len(seq) < percent:
				# Remove sequences that do not meet the threshold
				failed.append(species)
				with open(log, "a") as logfile:
					logfile.write(geneid + "\t" + species + 
								"\tLow Nucleotide Count\n")
		except ZeroDivisionError:
			failed.append(species)
	for i in failed:
		del newseq[i]
		n -= 1
	if n >= 2:
		return True, newseq, n
	elif n < 2:
		# Do not proceed if dict does not have at least two sequences
		pass 

#-----------------------------------------------------------------------------
	 
def writeDict(outfile, newseq):
	# Convert dictionary values to string and write to file
	with open(outfile, "w") as output:	 
		for species in newseq:
			output.write(">" + str(species) + "\n")
			seq = ""
			for codon in newseq[species]:
				seq += str(codon)
			output.write(seq + "\n")

def main():
	parser = argparse.ArgumentParser(description="This script will filter a \
multi-fasta file based on a specified reference species.")
	parser.add_argument("-i", help="Path to input directory.")
	parser.add_argument("-o", help="Path to output directory.")
	parser.add_argument("-r", help="Genome build name of the reference \
species as it appears in the fasta alignment.")
	parser.add_argument("-p", type=float, default=0.5,
help="Minimum required percentage of nucleotides remaining after filtering \
(as a decimal).")
	parser.add_argument("--retainStops", action="store_true", 
help="Specifies that sequences containing internal stop codons should be \
ratained.")
	# Parse arguments and assign to variables
	args = parser.parse_args()
	indir = args.i
	if indir[-1] != "/":
		indir += "/"
	outdir = args.o
	if outdir[-1] != "/":
		outdir += "/"
	ref = args.r
	percent = args.p
	retainStops = args.retainStops
	openFiles(indir, outdir, ref, percent, retainStops)

if __name__ == "__main__":
	main()
