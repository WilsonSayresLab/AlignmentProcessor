'''This program will run CodeML on a directory of single gene alignments.
It will generate a unique control file and tree file for each input gene
before invoking CodeML using the number of CPUs specified by the user
(default =1).

	Copyright 2016 by Shawn Rupp'''

from datetime import datetime
from sys import argv
from glob import glob
from subprocess import Popen
from shlex import split
from functools import partial
from multiprocessing import cpu_count
from multiprocessing.dummy import Pool
import shutil
import os

# Define max number of threads
MAXCPU = cpu_count()

def controlFiles(path):
	'''Reads input files and stores them in memory'''
	usertree = False
	# Make temp directory
	try:
		os.mkdir(path + "tmp/")
	except FileExistsError:
		pass
	with open(path + "codeml.ctl", "r") as control:
		ctl = control.readlines()
		for line in ctl:
			# Determine if a phylogenic tree is needed
			if "runmode = 0" in line or "runmode = 1" in line:
				usertree = True
	# Only run if a tree is required
	if usertree == True:
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
		with open(path + "tmp/ref.tree", "w") as reftree:
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
	elif usertree == False:
		pass
	return usertree, nodes, ctl, species

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
	'''Identifies genes which have already been run through CodeML.'''
	finished = path + "Logs/finishedCodeML.txt"
	if os.path.isfile(finished) == False:
		with open(finished, "w") as fin:
			# Create log file and empty list
			completed = []
	elif os.path.isfile(finished) == True:
		# Create list of completed files
		with open(finished, "r") as fin:
			completed = fin.readlines()
	return completed

def runCodeml(ap, usertree, path, completed, ctl, nodes, species, genes, gene):
	'''Creates temporary control and tree files and runs CodeML.'''
	go = False
	try:
		geneid = gene.split("/")[-1].split(".")[0]
		if (geneid + "\n") in completed:
			pass
		else:
			rmseqs = removedSeq(species, gene)
			# Create unique working directory (allows CodeML to create
			# multiple temp files with the same name).
			wd = path + "tmp/" + geneid
			try:
				os.mkdir(wd)
			except FileExistsError:
				pass
			# Set unique file names
			outfile = path + "07_codeml/" + geneid + ".mlc"
			tempctl = wd  + "/" + geneid + ".ctl"
			tmpTree = wd  + "/" + geneid + ".tree"
			if usertree == True:
				# Prune tree only if it will be used
				pruneTree(ap, path, geneid, tmpTree, rmseqs, nodes)
			with open(tempctl, "w") as temp:
				# creates unique control file
				for line in ctl:
					if "seqfile" in line:
						temp.write("\tseqfile = " + gene + "\n")
					elif "outfile" in line:
						temp.write("\toutfile = " + outfile + "\n")
					elif "treefile" in line:
						if usertree == True:
							temp.write("\ttreefile = " + tmpTree + "\n")
					else:
						temp.write(line)
			# Calls CodeML
			os.chdir(wd)
			f = open(os.devnull, "w")
			cm = Popen(split(ap + "/paml/bin/codeml " + tempctl), stdout = f)
			cm.wait()
			# Delete temp files and add to count when finished
			shutil.rmtree(wd)
			# Append gene ID to list of finishedCodeML.txt
			with open(path + "Logs/finishedCodeML.txt", "a") as finished:
				finished.write(geneid + "\n")
			# Print progress
			n = genes.index(gene)
			if n%2 == 0:
				print("\tProcessed", str(n), "of",  str(len(genes)), 
					"files.\n")
	except TypeError:
		# Skip entires with NoneType
		pass

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

def pruneTree(ap, path, geneid, tmpTree, rmseqs, nodes):
	'''Calls the ape package in R to remove species whose sequences have been
removed due to low ccontent.'''
	# Determine if file has removed sequences
	if rmseqs:
		# Call ape package to remove species from tree
		outtree = path + "07_codeml/" + geneid + ".pruned"
		exclude = ""
		for i in rmseqs:
			exclude += i
		# Call R script
		pt = Popen(split("Rscript " + ap + "bin/07_pruneTree.R " +
					path + "tmp/ref.tree " + outtree + " " + exclude))
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

def main():
	if argv[1] == "-h" or argv[1] == "--help":
		print("Usage: python 07_CodeMLonDir.py <path to CodeML control file> \
- <# of CPUs>")
		quit()
	else:
		starttime = datetime.now()
		# Save path to the AlignmentProcessor directory
		ap = os.getcwd() + "/"
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
		usertree, nodes, ctl, species = controlFiles(path)
		completed = outputFiles(path)
		# Call CodeML for all files in a directory.
		genes = glob(path + "06_phylipFiles/" + "*.phylip")
		pool = Pool(processes = cpu)
		func = partial(runCodeml, ap, usertree, path, completed, ctl, 
						nodes, species, genes)
		print("\n\tRunning CodeML with", str(cpu), "threads....\n")
		try:
			pool.map(func, genes)
		except FileNotFoundError:
			pass		
		pool.close()
		pool.join()
		# Remove tmp directory
		print("\tCodeML runtime: ", datetime.now() - starttime, "\n")
		os.remove(path + "tmp/ref.tree")
		os.rmdir(path + "tmp/")

if __name__ == "__main__":
	main()
