'''This program will run CodeML on a directory of single gene alignments.
It will generate a unique control file and tree file for each input gene
before invoking CodeML using the number of CPUs specified by the user
(default =1).

	Copyright 2016 by Shawn Rupp'''

from datetime import datetime
from sys import stdout
from glob import glob
from subprocess import Popen
from shlex import split
from functools import partial
from multiprocessing import Pool, cpu_count
import argparse
import shutil
import os
import re

# Define max number of threads and devnull for capturing stdout
MAXCPU = cpu_count()
DEVNULL = open(os.devnull, "w")

#-----------------------------------------------------------------------------

def outputFiles(outdir):
	'''Identifies genes which have already been run through CodeML.'''
	finished = ""
	path = outdir.split("/")[:-2]
	for i in path:
		finished += i + "/"
	finished += "finishedCodeML.txt"
	if os.path.isfile(finished) == False:
		with open(finished, "w") as fin:
			# Create log file and empty list
			completed = []
	elif os.path.isfile(finished) == True:
		# Create list of completed files
		with open(finished, "r") as fin:
			completed = fin.readlines()
	return finished, completed

def controlFiles(indir, outdir, forward, completed):
	'''Reads input files and stores them in memory'''
	go = False
	pairwise = False
	multiple = False
	# Make temp directory
	try:
		os.mkdir(outdir + "tmp/")
	except FileExistsError:
		pass
	path = outdir.split("/")[:-2]
	out = ""
	for i in path:
		out += i + "/"
	control = glob(out + "*.ctl")
	if len(control) > 1:
		# Quit if multiple .ctl files are present
		print("\n\tPlease provide only one control file for CodeML.\n")
		quit()
	with open(control[0], "r") as infile:
		ctl = infile.readlines()
	for line in ctl:
		# Determine if a phylogenic tree is needed
		if "runmode = 0" in line or "runmode = 1" in line:
			go = makeTree(indir, outdir, forward, completed, ctl)
			multiple = True
		elif "runmode = -2" in line:
			go = pairwiseControl(indir, outdir, completed, ctl)
			pairwise = True
	if go == True:
		return pairwise, multiple

#-----------------------------------------------------------------------------

def makeTree(indir, outdir, forward, completed, ctl):
	'''Calls PhyML to create a gene tree.'''
	print("\tRunnning PhyMl to create gene trees...")
	genes = glob(indir + "*.phylip")
	for gene in genes:
		filename = gene.split("/")[-1]
		geneid = filename.split(".")[0]
		if (geneid + "\n") in completed:
			# Skip genes which have finished CodeML
			pass
		elif filename.split(".")[1] == "2":
			# Skip pairwise genes (PhyMl will not make a tree)
			pass
		else:
			# Create temp directory
			wd = outdir + "tmp/" + geneid + "/"
			try:
				os.mkdir(wd)
			except FileExistsError:
				pass
			# Set unique file names
			outfile = (outdir + geneid +"."+ filename.split(".")[1] + ".mlc")
			tempctl = wd + geneid + ".ctl"
			treefile = wd + filename + "_phyml_tree.txt"
			# Make unique control file
			makeCtl(gene, outfile, tempctl, treefile, ctl)
			# Call PhyML to make gene tree
			phy = Popen(split("PhyML/PhyML -q -i " + gene), stdout = DEVNULL)
			phy.wait()
			if phy.returncode == 0:
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
							# Find end of branch length
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

def pairwiseControl(indir, outdir, completed, ctl):
	'''Makes control files for pairwsie comparisons.'''
	print("\tMaking pairwise CodeML control files...")
	genes = glob(indir + "*.phylip")
	for gene in genes:
		filename = gene.split("/")[-1]
		geneid = filename.split(".")[0]
		if (geneid + "\n") in completed:
			# Skip genes which have finished CodeML
			pass
		else:
			# Create temp directory
			wd = outdir + "tmp/" + geneid + "/"
			try:
				os.mkdir(wd)
			except FileExistsError:
				pass
			# Set unique file names
			outfile = (outdir + geneid +"."+ filename.split(".")[1] + ".mlc")
			tempctl = wd + geneid + ".ctl"
			treefile = ""
			# Make unique control file
			makeCtl(gene, outfile, tempctl, treefile, ctl)
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

def runCodeml(ap, outdir, finished, completed, gene):
	'''Creates temporary control and tree files and runs CodeML.'''
	filename = gene.split("/")[-1]
	geneid = filename.split(".")[0]
	wd = outdir + "tmp/" + geneid + "/"
	if len(glob(wd + "*")) > 1:
		if filename.split(".")[1] == "2":
			# Skip pairwise genes if tree files are present
			pass
	if (geneid + "\n") in completed:
		pass
	else:
		tempctl = wd + geneid + ".ctl"
		# Calls CodeML
		os.chdir(wd)
		cm = Popen(split(ap + "/paml/bin/codeml " + tempctl), 
						stdout = DEVNULL)
		cm.wait()
		# Append gene ID to list of finishedCodeML.txt
		with open(finished, "a") as fin:
			fin.write(geneid + "\n")

#-----------------------------------------------------------------------------

def main():
	starttime = datetime.now()
	# Save path to the AlignmentProcessor directory
	ap = os.getcwd() + "/"
	run = False
	# Parse command
	parser = argparse.ArgumentParser(description="Runs CodeML on all files \
in a directory.")
	parser.add_argument("-i", help="Path to input directory.")
	parser.add_argument("-o", help="Path to output directory.")
	parser.add_argument("-t", type=int, default=1, help="Number of threads.")
	parser.add_argument("-f", default="", 
help="Forward species (name must be the same as it appears in input files.")
	parser.add_argument("--noCleanUp", action="store_false", 
help="Keep temporary files.")
	args = parser.parse_args()
	# Assign arguments
	indir = args.i
	if indir[-1] != "/":
		indir += "/"
	outdir = args.o
	if outdir[-1] != "/":
		outdir += "/"
	cpu = args.t
	if cpu > MAXCPU:
		cpu = MAXCPU
	forward = args.f
	cleanup = args.noCleanUp
	# Reads in required data
	finished, completed = outputFiles(outdir)
	pairwise, multiple = controlFiles(indir, outdir, forward, completed)
	if pairwise == True or multiple == True:
		# Call CodeML after PhyML completes.
		genes = glob(indir + "*.phylip")
		l = int(len(genes))
		pool = Pool(processes = cpu)
		func = partial(runCodeml, ap, outdir, finished, completed)
		print("\tRunning CodeML with", str(cpu), "threads....")
		rcml = pool.imap(func, genes, chunksize = int(l/cpu))
		pool.close()
		pool.join()	
	# Remove tmp directory
	if cleanup == True:
		shutil.rmtree(outdir + "tmp/")
	print("\tCodeML runtime: ", datetime.now() - starttime, "\n")

if __name__ == "__main__":
	main()
