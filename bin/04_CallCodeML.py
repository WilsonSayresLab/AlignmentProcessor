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
	cmd = ("python bin/04_callPhyML.py -i" + indir + " -o " + tmp + 
			" -t " + str(cpu))
	if forward:
		cmd += " -f " + forward
	phyml = Popen(split(cmd))
	phyml.wait()
	#if phyml.returncode() == 0:
	return True

#-----------------------------------------------------------------------------

def runCodeml(ap, outdir, finished, completed, gene):
	'''Creates temporary control and tree files and runs CodeML.'''
	filename = gene.split("/")[-1]
	geneid = filename.split(".")[0]
	wd = outdir + "tmp/" + geneid + "/"
	if filename.split(".")[1] == "2":
		if len(glob(wd + "*")) > 1:
			# Skip pairwise genes if tree files are present
			pass
	elif (geneid + "\n") in completed:
		pass
	else:
		tempctl = wd + geneid + ".ctl"
		# Calls CodeML
		os.chdir(wd)
		cm = Popen(split(ap + "/paml/bin/codeml " + tempctl), 
						stdout = DEVNULL)
		cm.wait()
		if cm.returncode == 0:
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
	run = controlFiles(indir, outdir, forward, cpu)
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
		func = partial(runCodeml, ap, outdir, finished, completed)
		print("\tRunning CodeML with", str(cpu), "threads....")
		rcml = pool.imap(func, genes, chunksize = chunk)
		pool.close()
		pool.join()	
	# Remove tmp directory
	if cleanup == True:
		shutil.rmtree(outdir + "tmp/")
	print("\tCodeML runtime: ", datetime.now() - starttime, "\n")

if __name__ == "__main__":
	main()
