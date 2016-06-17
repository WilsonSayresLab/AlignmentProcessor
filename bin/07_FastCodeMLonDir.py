'''This program will run fastCodeML on a directory of single gene alignments.
It will generate a unique control file and tree file for each input gene
before invoking fastCodeML using the number of CPUs specified by the user
(default =1).

	Copyright 2016 by Shawn Rupp'''

from datetime import datetime
from sys import argv
from glob import glob
from subprocess import Popen
from shlex import split
from multiprocessing import cpu_count
import shutil
import os

# Define max number of threads
MAXCPU = cpu_count()

def controlFiles(path):
	'''Reads input files and stores them in memory'''
	with open(path + "codeml.tree", "r") as intree:
		tree = intree.readlines()[0].rstrip()
	# Create list of species and coresponding nodes
	nodes = {}
	splt = tree.split(",")
	for i in splt:
		i = i.replace("(", "")
		i = i.replace(")", "")
		i = i .replace(";", "")
		if "$" in i:
			print()
			print("\tPlease only use the pound sign (#) to indicate nodes.")
			print()
			quit()
		if "#" in i:
			node = i.split()
			nodes[node[0]] = node[1]
	with open(path + "07_codeml/ref.tree", "w") as reftree:
		# Add semicolon to the end of the tree if it is not present
		if tree[-1] != ";":
			tree += ";"
		# Remove node markers and write corrected tree to file
		if "#" in tree:
			for j in range(0, tree.count("#")):
				i = tree.index("#")
				tree = tree[:i-1] + tree[i+2:]
		reftree.write(tree)
	species = allSpecies(tree)
	return nodes, species

def allSpecies(tree):
	'''Determines all species in codeml.tree'''
	species = []
	for i in tree.split(","):
		# Isolate species name 
		j = i.strip().split()[0]
		j = j.replace("(", "")
		j = j.replace(")", "")
		j = j.replace(";", "")
		j = j.split(":")[0]
		species.append(j)
	return species

def outputFiles(path):
	'''Identifies genes which have already been run through fastCodeML.'''
	finished = path + "Logs/finishedFastCodeML.txt"
	if os.path.isfile(finished) == False:
		with open(finished, "w") as fin:
			# Create log file and empty list
			completed = []
	elif os.path.isfile(finished) == True:
		# Create list of completed files
		with open(finished, "r") as fin:
			completed = fin.readlines()
	return completed

def	runFCML(path, cpu, completed, nodes, species, genes):
	'''Creates temporary control and tree files and runs fastCodeML.'''
	for gene in genes:
		geneid = gene.split("/")[-1].split(".")[0]
		if (geneid + "\n") in completed:
			pass
		else:
			go = False
			rmseqs = removedSeq(species, gene)
			# Set unique file names
			outfile = path + "07_codeml/" + geneid + ".fcml"
			tmpTree = path + "07_codeml/" + geneid + ".tree"
			go = pruneTree(path, geneid, tmpTree, rmseqs, nodes)
			if go == True:
				# Calls fastCodeML
				f = open(os.devnull, "w")
				fcml = Popen(split("./FastCodeML/fast -bf -ps -nt " + str(cpu)
										 + " -ou " + outfile + " " + tmpTree + " "
										 + gene), stdout = f)
				fcml.wait()
				# Delete temp files and add to count when finished
				os.remove(tmpTree)
				# Append gene ID to list of finishedFastCodeML.txt
				with open(path + "Logs/finishedFastCodeML.txt", "a") as finished:
					finished.write(geneid + "\n")
				# Print progress
				n = genes.index(gene)
				if n%1000 == 0:
					print("\tProcessed", str(n), "of",  str(len(genes)), 
						"files.\n")
	return True

def removedSeq(species, gene):
	'''Identifies species which are missing from alignment'''
	# Determine which species remain in the alignment
	rmseqs = []
	remaining = []
	with open(gene, "r") as alignment:
		for line in alignment:
			if line[0] == " ":
				# Skip first line
				pass
			else:
				remaining.append(line.split()[0])
	for i in species:
		if i not in remaining:
			rmseqs.append(i)
	return rmseqs

def pruneTree(path, geneid, tmpTree, rmseqs, nodes):
	'''Calls the ape package in R to remove species whose sequences have been
removed due to low ccontent.'''
	# Determine if file has removed sequences
	if rmseqs:
		# Call ape package to remove species from tree
		first = True
		outtree = path + "07_codeml/" + geneid + ".pruned"
		for i in rmseqs:
			# Call R script
			if first == True:
				pt = Popen(split("Rscript bin/07_pruneTree.R " +
							path + "07_codeml/ref.tree " + outtree + " " + i))
				pt.wait()
				first = False
			elif first == False:
				pt = Popen(split("Rscript bin/07_pruneTree.R " +
							outtree + " " + outtree + " "  + i))
				pt.wait()
		# Read in ape output
		with open(outtree, "r") as pruned:
			tree = pruned.readlines()
		os.remove(outtree)
		# Add nodes back into tree
		for species in nodes:
			if species in tree:
				# Determine location and length of species name
				i = tree.index(species) + len(species)
				# Insert space and  node symbol after species name
				tree = (tree[:i] + " " + str(nodes[species]) +
						tree[i:])
		with open(tmpTree, "w") as temptree:
			# Write tree to file
			string = ""
			for i in tree:
				string += i
			temptree.write(string)			
	elif not rmseqs:
		# Read in codeml.tree instead of pruned tree
		with open(path + "codeml.tree", "r") as ref:
			tree = ref.readlines()[0]
		with open(tmpTree, "w") as temptree:
			temptree.write(tree)
	return True

def main():
	if argv[1] == "-h" or argv[1] == "--help":
		print("Usage: python 07_CodeMLonDir.py <path to input directory> \
- <# of CPUs>")
		quit()
	else:
		starttime = datetime.now()
		cpu = 1
		retainStops = False
		# Parse command
		for i in argv:
			if i == "-i":
				path = argv[argv.index(i) + 1]
			elif i == "-t":
				cpu = int(argv[argv.index(i) + 1])
		# Set directory names and add a trailing "/" if necessary
		if path[-1] != "/":
			path += "/"
		# Make sure too many threads have not been specified
		if cpu > MAXCPU:
			cpu = MAXCPU
		# Reads in required data
		nodes, species = controlFiles(path)
		completed = outputFiles(path)
		# Call fastCodeML for all files in a directory.
		genes = glob(path + "06_phylipFiles/" + "*.phylip")
		print("\n\tRunning FastCodeML with", str(cpu), "threads....\n")
		go = runFCML(path, cpu, completed, nodes, species, genes)
		if go == True:
			os.remove(path + "07_codeml/ref.tree")
			print("\tFastCodeML runtime: ", datetime.now() - starttime, "\n")

if __name__ == "__main__":
	main()
