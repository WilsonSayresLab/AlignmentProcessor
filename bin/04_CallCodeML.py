'''This program will run CodeML on a directory of single gene alignments.
It will generate a unique control file and tree file for each input gene
before invoking CodeML using the number of CPUs specified by the user
(default = 1).

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
import math
import os
import re
from callPhyML import phyml

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

def controlFiles(indir, outdir, forward, cpu):
	'''Reads input files and stores them in memory'''
	# Make temp directory
	tmp = outdir + "tmp/"
	try:
		os.mkdir(tmp)
	except FileExistsError:
		pass
	multiple = phyml(indir, tmp, cpu, forward)
	return True, multiple

#-----------------------------------------------------------------------------

def runCodeml(ap, outdir, finished, completed, multiple, gene):
	'''Creates temporary control and tree files and runs CodeML.'''
	filename = gene.split("/")[-1]
	geneid = filename.split(".")[0]
	wd = outdir + "tmp/" + geneid + "/"
	if (geneid + "\n") in completed:
		pass
	else:
		tempctl = wd + "codeml.ctl"
		os.chdir(wd)
		if multiple == True:
			if filename.split(".")[1] == "2":
				pass
			else:
				# Calls CodeML if 3 or more sequences are present
				cm = Popen(split(ap + "paml/bin/codeml " + tempctl), 
								stdout = DEVNULL)
		elif multiple == False:
			# Call CodeML for all files
			cm = Popen(split(ap + "paml/bin/codeml " + tempctl), 
							stdout = DEVNULL)
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
	parser.add_argument("--cleanUp", action="store_true", 
help="Remove temporary files (it may be useful to retain phylogenic trees for future use).")
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
	cleanup = args.cleanUp
	# Reads in required data
	finished, completed = outputFiles(outdir)
	run, multiple = controlFiles(indir, outdir, forward, cpu)
	if run == True:
		# Call CodeML after PhyML completes.
		genes = glob(indir + "*.phylip")
		l = int(len(genes))
		# Determine chunksize
		if l <= cpu:
			chunk = 1
		elif l > cpu:
			chunk = int(math.ceil(l/cpu))
		pool = Pool(processes = cpu)
		func = partial(runCodeml, ap, outdir, finished, completed, multiple)
		print(("\tRunning CodeML with {0!s} threads....").format(cpu))
		rcml = pool.imap_unordered(func, genes, chunksize = chunk)
		pool.close()
		pool.join()		
	# Remove tmp directory
	if cleanup == True:
		shutil.rmtree(outdir + "tmp/")
	print(("\tCodeML runtime: {0!s}").format(datetime.now() - starttime))

if __name__ == "__main__":
	main()
