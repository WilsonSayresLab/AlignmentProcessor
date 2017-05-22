'''This program will run CodeML on a directory of single gene alignments.
It will generate a unique control file and tree file for each input gene
before invoking CodeML using the number of CPUs specified by the user
(default = 1).

	Copyright 2016 by Shawn Rupp'''

from sys import stdout
from glob import glob
from subprocess import Popen
from shlex import split
from shutil import rmtree
import math
import shutil
import os
import re

DEVNULL = open(os.devnull, "w")

def makeTree(ap, gene, wd, treefile, forward):
	'''Calls PhyML to create a gene tree.'''
	# Call PhyML to make gene tree
	phy = Popen(split(ap + "PhyML/PhyML -q -i " + gene), stdout = DEVNULL)
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

def parallelize(ap, outdir, finished, completed, multiple, cpu, ctl, 
				forward, treefile, gene):
	filename = gene.split("/")[-1]
	geneid = filename.split(".")[0]
	if (geneid + "\n") in completed:
		# Skip genes which have been analyzed
		pass
	else:
		# Make working directory and gene control file
		wd = outdir + "tmp/" + geneid + "/"
		try:
			os.mkdir(wd)
		except FileExistsError:
			pass
		tempctl = wd + "codeml.ctl"
		outfile = (outdir + geneid +"."+ filename.split(".")[1] + ".mlc")
		if multiple == True:
			if filename.split(".")[1] == "2":
				# Skip pairwise genes and remove directory
				rmtree(wd)
				pass
			else:
				if not treefile:
					# Run Phyml
					treefile = wd + filename + "_phyml_tree.txt"
					makeTree(ap, gene, wd, treefile, forward)
				# Make control file
				makeCtl(gene, outfile, tempctl, treefile, ctl)
				os.chdir(wd)
				# Call CodeML
				with open(wd + "codemlLog.txt", "w") as tmpout:
					cm = Popen(split(ap + "paml/bin/codeml " + tempctl),
									shell = True, stdout = tmpout)
					cm.wait()
		elif multiple == False:
			# Make control file
			treefile = wd + filename + ""
			makeCtl(gene, outfile, tempctl, treefile, ctl)
			# Call CodeML for all files
			os.chdir(wd)
			with open(wd + "codemlLog.txt", "w") as tmpout:
				cm = Popen(split(ap + "paml/bin/codeml " + tempctl),
								shell = True, stdout = tmpout)
				cm.wait()
		with open(finished, "a") as fin:
			# Append gene id to list when done
			fin.write(geneid + "\n")
