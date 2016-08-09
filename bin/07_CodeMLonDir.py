'''This program will run CodeML on a directory of single gene alignments.
It will generate a unique control file and tree file for each input gene
before invoking CodeML using the number of CPUs specified by the user
(default =1).

	Copyright 2016 by Shawn Rupp'''

from datetime import datetime
from sys import argv, stdout
from glob import glob
from subprocess import Popen
from shlex import split
from functools import partial
from multiprocessing import Pool, cpu_count
import shutil
import os
import re

# Define max number of threads and devnull for capturing stdout
MAXCPU = cpu_count()
DEVNULL = open(os.devnull, "w")

#-----------------------------------------------------------------------------

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

def controlFiles(path, forward, completed):
	'''Reads input files and stores them in memory'''
	go = False
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
			go = makeTree(path, forward, completed, ctl)
	if go == True:
		return ctl

def makeTree(path, forward, completed, ctl):
	'''Calls PhyML to create a gene tree.'''
	print("\n\tRunnning PhyMl to create gene trees...\n")
	genes = glob(path + "06_phylipFiles/" + "*.phylip")
	for gene in genes:
		filename = gene.split("/")[-1]
		geneid = filename.split(".")[0]
		if (geneid + "\n") in completed:
			# Skip genes which have finished CodeML
			pass
		else:
			# Create temp directory
			wd = path + "tmp/" + geneid + "/"
			try:
				os.mkdir(wd)
			except FileExistsError:
				pass
			# Set unique file names
			outfile = (path + "07_codeml/" + geneid + "." +
						filename.split(".")[1] + ".mlc")
			tempctl = wd + geneid + ".ctl"
			treefile = wd + filename + "_phyml_tree.txt"
			# Make unique control file
			makeCtl(gene, outfile, tempctl, treefile, ctl)
			# Call PhyML to make gene tree
			phy = Popen(split("PhyML/PhyML -q -i " + gene), stdout = DEVNULL)
			phy.wait()
			# Move PhyML output to temp directory
			output = glob(gene + "_phyml_*") 
			for i in output:
				shutil.copy(i, wd)
				os.remove(i)
			# Read in PhyML tree
			with open(treefile, "r") as genetree:
				try:
					tree = genetree.readlines()[0]
				except IndexError:
					pass
			# Remove branch lables introduced by PhyML
			tree = re.sub(r"\d+\.\d+:", ":", tree)
			# Add forward node to tree if specified 
			if forward:
				if forward in tree:
					# Determine location and length of species name
					i = tree.index(forward) + len(forward)
					if ":" in tree:
						# Find end of branch length (find next end-of-secies char)
						comma = tree.find(",", i)
						paren = tree.find(")", i)
						i = min([comma, paren])
					# Insert space and node symbol after species name
					tree = (tree[:i] + " #1" + tree[i:])
			elif forward not in tree:
				pass
			with open(treefile, "w") as outtree:
				# Overwrite treefile
				string = ""
				for i in tree:
					string += i
				outtree.write(string)	
	return True		

def makeCtl(gene, outfile, tempctl, treefile, ctl):
	'''Creates unique control file'''
	with open(tempctl, "w") as temp:
		for line in ctl:
			if "seqfile" in line:
				temp.write("\tseqfile = " + gene + "\n")
			elif "outfile" in line:
				temp.write("\toutfile = " + outfile + "\n")
			elif "treefile" in line:
				temp.write("\ttreefile = " + treefile + "\n")
			else:
				temp.write(line)

#-----------------------------------------------------------------------------

def runCodeml(ap, path, completed, gene):
	'''Creates temporary control and tree files and runs CodeML.'''
	filename = gene.split("/")[-1]
	geneid = filename.split(".")[0]
	if (geneid + "\n") in completed:
		pass
	else:
		# Create unique working directory (allows CodeML to create
		# multiple temp files with the same name).
		wd = path + "tmp/" + geneid + "/"
		tempctl = wd + geneid + ".ctl"
		# Calls CodeML
		os.chdir(wd)
		cm = Popen(split
			(ap + "/paml/bin/codeml " + tempctl), stdout = DEVNULL)
		cm.wait()
		# Delete temp files and add to count when finished
		shutil.rmtree(wd)
		# Append gene ID to list of finishedCodeML.txt
		with open(path + "Logs/finishedCodeML.txt", "a") as finished:
			finished.write(geneid + "\n")

#-----------------------------------------------------------------------------

def main():
	if argv[1] == "-h" or argv[1] == "--help":
		print("Usage: python 07_CodeMLonDir.py -i \
<path to CodeML control file> -t <# of CPUs> -f forward species node")
		quit()
	else:
		starttime = datetime.now()
		# Save path to the AlignmentProcessor directory
		ap = os.getcwd() + "/"
		cpu = 1
		retainStops = False
		forward = ""
		# Parse command
		for i in argv:
			if i == "-i":
				path = argv[argv.index(i) + 1]
			elif i == "-t":
				cpu = int(argv[argv.index(i) + 1])
			elif i == "-f":
				forward = argv[argv.index(i) + 1]
		# Set directory names and add a trailing "/" if necessary
		if path[-1] != "/":
			path += "/"
		# Make sure too many threads have not been specified
		if cpu > MAXCPU:
			cpu = MAXCPU
		# Reads in required data
		completed = outputFiles(path)
		ctl = controlFiles(path, forward, completed)
		if ctl:
			# Call CodeML after PhyML completes.
			genes = glob(path + "06_phylipFiles/" + "*.phylip")
			l = int(len(genes))
			pool = Pool(processes = cpu)
			func = partial(runCodeml, ap, path, completed)
			print("\n\tRunning CodeML with", str(cpu), "threads....\n")
			rcml = pool.imap(func, genes, chunksize = int(l/cpu))
			pool.close()
			pool.join()	
			# Remove tmp directory
			os.rmdir(path + "tmp/")
			print("\n\tCodeML runtime: ", datetime.now() - starttime, "\n")

if __name__ == "__main__":
	main()
