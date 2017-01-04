'''This script will run PhyML on a directory of single gene alignments.
It will generate a unique control file and tree file for each input gene.

	Copyright 2016 by Shawn Rupp'''

from sys import stdout
from glob import glob
from subprocess import Popen
from shlex import split
from functools import partial
from multiprocessing import Pool, cpu_count
import argparse
import shutil
import math
import os
import re

DEVNULL = open(os.devnull, "w")

#-----------------------------------------------------------------------------

def outputFiles(outdir):
	'''Identifies genes which have already been run through PhyML.'''
	completed = []
	multiple = False
	tmpdirs = glob(outdir + "*")
	for i in tmpdirs:
		# Extract gene IDs and append to list
		completed.append(i.split("/")[-1].replace("/", ""))
	# Reconstruct output path
	path = outdir.split("/")[:-3]
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
			multiple = True
	return completed, ctl, multiple

#-----------------------------------------------------------------------------

def makeTree(indir, outdir, forward, completed, ctl, gene):
	'''Calls PhyML to create a gene tree.'''
	filename = gene.split("/")[-1]
	geneid = filename.split(".")[0]
	if filename.split(".")[1] == "2":
		# Skip pairwise genes (PhyMl will not make a tree)
		pass
	else:
		# Create temp directory and codeml ouput directory
		path = outdir.split("/")[:-2]
		out = ""
		for i in path:
			out += i + "/"
		wd = outdir + geneid + "/"
		try:
			os.mkdir(wd)
		except FileExistsError:
			pass
		# Set unique file names
		outfile = (out + geneid + "." + filename.split(".")[1] + ".mlc")
		tempctl = wd + "codeml.ctl"
		treefile = wd + filename + "_phyml_tree.txt"
		# Make unique control file
		makeCtl(gene, outfile, tempctl, treefile, ctl)
		if (geneid) in completed:
			# Skip genes which have finished CodeML
			pass
		else:
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

def pairwiseControl(indir, outdir, ctl, gene):
	'''Makes control files for pairwsie comparisons.'''
	filename = gene.split("/")[-1]
	geneid = filename.split(".")[0]
	# Create temp directory
	wd = outdir + geneid + "/"
	try:
		os.mkdir(wd)
	except FileExistsError:
		pass
	# Create temp directory
	path = outdir.split("/")[:-2]
	out = ""
	for i in path:
		out += i + "/"
	# Set unique file names
	outfile = (out + geneid +"."+ filename.split(".")[1] + ".mlc")
	tempctl = wd + "codeml.ctl"
	treefile = ""
	# Make unique control file
	makeCtl(gene, outfile, tempctl, treefile, ctl)

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

def phyml(indir, outdir, cpu, forward=""):
	# Reads in required data
	completed, ctl, multiple = outputFiles(outdir)
	# Call PhyML for multiple alignments or write pairwise control files.
	genes = glob(indir + "*.phylip")
	l = int(len(genes))
	# Determine chunksize
	if l <= cpu:
		chunk = 1
	elif l > cpu:
		chunk = int(math.ceil(l/cpu))
	pool = Pool(processes = cpu)
	if multiple == True:
		# Call PhyML and make multiple control files
		print(("\tRunnning PhyMl with {0!s} threads to create gene \
trees...").format(cpu))
		func = partial(makeTree, indir, outdir, forward, completed, ctl)
		rpml = pool.imap(func, genes, chunksize = chunk)
		pool.close()
		pool.join()	
	elif multiple == False:
		# Make pairwise control files
		print(("\tMaking pairwise CodeML control files with {0!s} \
threads...").format(cpu))
		func = partial(pairwiseControl, indir, outdir, ctl)
		rpml = pool.imap(func, genes, chunksize = chunk)
		pool.close()
		pool.join()	
	return multiple
