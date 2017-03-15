'''This program will call the CodeML cython script on a directory of 
single gene alignments. It will generate a unique control file and tree 
file for each input gene before invoking CodeML using the number of CPUs 
specified by the user (default = 1).

	Copyright 2016 by Shawn Rupp'''

from __future__ import division
from datetime import datetime
from glob import glob
from multiprocessing import Pool, cpu_count
from functools import partial
import argparse
import shutil
import sys
import math
import os
from parallelCodeML import parallelize

MAXCPU = cpu_count()

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
	multiple = False
	# Make temp directory
	tmp = outdir + "tmp/"
	try:
		os.mkdir(tmp)
	except FileExistsError:
		pass
		# Reconstruct output path
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
			multiple = True
	return ctl, multiple

#-----------------------------------------------------------------------------

def main():
	starttime = datetime.now()
	# Save path to the AlignmentProcessor directory
	ap = os.getcwd() + "/"
	if " " in ap:
		# Change to warning ########################################################
		print("\tWARNING: AlignmentProcessor will not run properly if there \
is a space in its PATH name.")
		ap = ap.replace(" (ASU)", "")
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
help="Remove temporary files (it may be useful to retain phylogenic trees \
for future use).")
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
	# Reads in required data
	finished, completed = outputFiles(outdir)
	ctl, multiple = controlFiles(indir, outdir, forward, cpu)
	# Call PhyML and CodeML in parallel completes.
	if ctl:
		# Call CodeML and PhyML
		genes = glob(indir + "*.phylip")
		l = int(len(genes))
		func = partial(parallelize, ap, outdir, finished, completed, multiple,
						cpu, ctl, forward)
		print(("\tRunning CodeML on {0!s} genes with {1!s} threads...."
				).format(l, cpu))
		pool = Pool(processes = cpu)
		for i, _ in enumerate(pool.imap_unordered(func, genes), 1):
			sys.stderr.write("\r\t{0:%} of genes have finished".format(i/l))
		pool.close()
		pool.join()
	# Remove tmp directory
	if args.cleanUp == True:
		shutil.rmtree(outdir + "tmp/")
	print(("\n\tCodeML runtime: {0!s}").format(datetime.now() - starttime))

if __name__ == "__main__":
	main()
