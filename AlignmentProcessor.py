'''AlignmentProcessor will run the subsituion rate pipeline to produce trimmed
axt or phylip files for use with KaKs_calculator or PhyMl.

	Copyright 2016 by Shawn Rupp

	This package is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation version 3 of the License.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.'''

from datetime import datetime
import argparse
from subprocess import Popen
from shlex import split
from glob import glob
import os

def checkInput(ref, kaks, codeml, outdir):
	'''Makes sure necessary programs are installed and proper file conversion
is run for optional analysis program.'''
	proceed = True
	if not ref:
		print("\n\tError: Please specify a reference species.\n")
		proceed = False
	if kaks == True:
		if codeml == True:
			print("\n\tError: Please specify only one analysis program.\n")
			proceed = False
		indir = os.path.isfile("bin/KaKs_Calculator")
		if indir == False:
			print("\n\tError: Please install KaKs_Calculator in the\
 AlignmentProcessor bin.\n")
			proceed = False
	if codeml == True:
		control = glob(outdir + "*.ctl")
		if len(control) == 0:
			print("\n\tPlease supply a control file for CodeML.\n")
			proceed = False
		elif len(control) > 1:
			print("\n\tPlease supply only one CodeML control file.\n")
			proceed = False
		indir = os.path.isfile("paml/bin/codeml")
		phyml = os.path.isfile("PhyML/PhyML")
		if indir == False:
			print("\n\tError: Please install PAML in the\
 AlignmentProcessor folder.\n")
			proceed = False
		if phyml == False:
			print("\n\tError: Please install and rename PhyML for use\
 with CodeML.\n")
			proceed = False
	return proceed

		
def makeDir(path, outdir, axt, phylip, kaks, codeml):
	'''Makes all sub-directories used by program.'''
	print("\n\tMaking output directories...")
	os.chdir(outdir)
	for i in ["01_splitFasta", "02_filteredFasta"]:
		try:
			os.mkdir(i)
		except FileExistsError:
			pass
	if axt == True:
		try:
			os.mkdir("03_axtFiles")
		except FileExistsError:
			pass
	if kaks == True:
		try:
			os.mkdir("04_KaKsOutput")
		except FileExistsError:
			pass
	if phylip == True:
		try:
			os.mkdir("03_phylipFiles")
		except FileExistsError:
			pass
	if codeml == True:
		try:
			os.mkdir("04_CodemlOutput")
		except FileExistsError:
			pass	
	os.chdir(path)

#-----------------------------------------------------------------------------
			
def splitFasta(infile, outdir):
	'''Splits fasta alignment into one file per gene.'''
	sf = Popen(split("python bin/01_SplitFasta.py -i "
				+ infile + " -o " + outdir + "01_splitFasta"))
	sf.wait()
	if sf.returncode == 0:
		return True
	
def filterFasta(outdir, ref, percent, retainstops):
	'''Calls fasta filtering script'''
	cmd = ("python bin/02_FilterFasta.py -i " + outdir + "01_splitFasta -o "
			 + outdir + "02_filteredFasta -r " + ref)
	if retainstops == True:
		cmd += " --retainStops"
	ff = Popen(split(cmd))
	ff.wait()
	if ff.returncode == 0:
		return True
	
def convert(outdir, axt, phylip):
	'''Converts fasta files to specified output type'''
	cmd = ("python bin/03_ConvertFasta.py -i " + outdir + "02_filteredFasta -o "
			+ outdir)
	if phylip == True:
		cmd += "03_phylipFiles --phylip"
	elif axt == True:
		cmd += "03_axtFiles --axt"
	cf = Popen(split(cmd))
	cf.wait()
	if cf.returncode == 0:
		return True

def calculateKaKs(outdir, method, cpu):
	'''Calls KaKs_Calculator to calculate substition rates.'''
	ck = Popen(split("python bin/04_CallKaKs.py -i" + outdir + 
		"03_axtFiles -o " + outdir + "04_KaKsOutput -m " + method + 
		" -t  " + str(cpu)))
	ck.wait()
	if ck.returncode == 0:
		return True

def runcodeml(cpu, outdir, forward, ntree, cleanup):
	'''Runs codeml on a directory.'''
	# Build commands and add options if necessary
	cmd = ("python bin/04_CallCodeML.py -t " + str(cpu) + " -i " + outdir + 
			"03_phylipFiles" + " -o " + outdir + "04_CodemlOutput")
	if cleanup == True:
		cmd += " --cleanUp"
	if forward:
		cmd += " -f " + forward
	if ntree:
		cmd += " -n " + ntree
	cm = Popen(split(cmd))
	cm.wait()
	if cm.returncode == 0:
		return True

#-----------------------------------------------------------------------------

def main():
	starttime = datetime.now()
	# Set arguments
	parser = argparse.ArgumentParser(description="AlignmentProcessor will run \
the subsituion rate pipeline to produce trimmed axt or phylip files for use \
with KaKs_calculator or PhyMl.\nAlignmentProcessor1.6 Copyright 2016 by \
Shawn Rupp\nThis program comes with ABSOLUTELY NO WARRANTY\nThis is free \
software, and you are welcome to redistribute it under certain conditions\n")
	parser.add_argument("-i", help="Path to input file.")
	parser.add_argument("-o", help="Path to output directory.")
	parser.add_argument("-r", help="Genome build name of the reference \
species as it appears in the fasta alignment.")
	parser.add_argument("-p", type=float, default=0.5,
help="Minimum required percentage of nucleotides remaining after filtering \
(as a decimal).")
	parser.add_argument("--retainStops", action="store_true", 
help="Specifies that sequences containing internal stop codons should be \
ratained.")
	parser.add_argument("--axt", action="store_true", 
help="Convert files to axt format without calling KaKs_Calculator.")
	parser.add_argument("--phylip", action="store_true",
help="Convert files to sequential phylip format without calling CodeML.")
	parser.add_argument("--kaks", action="store_true",
help="Calls KaKs_Calculator on axt files.")
	parser.add_argument("-m", default="NG", 
help="Method for calculating Ka/Ks.")
	parser.add_argument("--codeml", action="store_true",
help="Calls CodeMl on phylip files.")
	parser.add_argument("-t", type=int, default=1, help="Number of threads.")
	parser.add_argument("-f", default="", 
help="Forward species (name must be the same as it appears in input files.")
	parser.add_argument("-n", help = "Path to optional Newick tree.")
	parser.add_argument("--cleanUp", action="store_true", 
help="Remove temporary files (it may be useful to retain phylogenic trees \
for future use).")
	# Parse arguments and assign to variables
	args = parser.parse_args()
	infile = args.i
	outdir = args.o
	if outdir != "/":
		outdir += "/"
	ref = args.r
	percent = args.p
	retainstops = args.retainStops
	axt = args.axt
	phylip = args.phylip
	kaks = args.kaks
	method = args.m
	codeml = args.codeml
	cpu = args.t
	forward = args.f
	ntree = args.n
	cleanup = args.cleanUp
	# Check input commands prior to running:
	proceed = checkInput(ref, kaks, codeml, outdir)
	if proceed == True:
		# Set axt or phylip
		if codeml == True:
			phylip = True
		elif kaks == True:
			axt = True
		# Set checkpoint variables to False:
		sf = False
		ff = False
		cf = False
		done = False
		# Save working directory to variable and call other scripts:
		path = os.getcwd()
		path = path + "/"
		makeDir(path, outdir, axt, phylip, kaks, codeml)
		sf = splitFasta(infile, outdir)
		if sf == True:
			ff = filterFasta(outdir, ref, percent, retainstops)
		if ff == True:
			cf = convert(outdir, axt, phylip)
		if cf == True:
			# Run KaKs_Calculator:
			if kaks == True:
				done = calculateKaKs(outdir, method, cpu)
			# Run codeml
			elif codeml == True:
				done = runcodeml(cpu, outdir, forward, ntree, cleanup)
			else:
				# Exit if neither program was called
				done = True
		# Print run time
		if done == True:
			print("\n\tTotal runtime: ", datetime.now() - starttime, "\n")

if __name__ == "__main__":
	main()
