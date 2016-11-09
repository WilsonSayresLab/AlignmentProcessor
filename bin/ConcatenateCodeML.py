'''This script will concatenate CodeML output from AlignmentProcessor
 and print them in a single file.'''

import argparse
from glob import glob

def concatenate(multiple, pairwise, indir, outfile):
	# Parses CodeML output and prints to file
	if multiple == True and pairwise == True:
		print("\n\tPlease specify only one alignment type.\n")
		quit()
	elif multiple == True:
		# Read in codeml output for multiple alignments
		readMultiple(indir, outfile)
	elif pairwise == True:
		# Read in codeml output for pairwise alignments
		readPairwise(indir, outfile)

def readMultiple(indir, outfile):
	with open(outfile, "w") as output:
		# Open output file and write header
		output.write("TranscriptID,dN,dS,dN/dS,TreeLength,lnL,Species\n")
		infiles = glob(indir + "*")
		for infile in infiles:
			# Get gene id and number of species remaining for each file
			filename = infile.split("/")[-1]
			transcript = filename.split(".")[0]
			species = filename.split(".")[1]
			if int(species) > 2:
				with open(infile, "r") as mlc:
					for line in mlc:
						# Get substitution rates, tree lengths, etc
						if "tree length =" in line:
							length = line.split("=")[1].strip()
						elif "tree length for dN" in line:
							dn = line.split(":")[1].strip()
						elif "tree length for dS" in line:
							ds = line.split(":")[1].strip()
						elif "lnL" in line:
							lnl = line.split("):")[1]
							lnl = lnl.split()[0].strip()
				# Calculate dN/dS and save as a string
				try:
					dnds = str(float(dn)/float(ds))
				except ZeroDivisionError:
					dnds = "NA"
				# Append data to list and convert into string
				data = [dn, ds, dnds, length, lnl, species]
				transcript += ","
				for i in data:
					transcript += str(i) + ","
				transcript = transcript[:-1]
				output.write(transcript + "\n")
			else:
				# Skip files with only two sequences
				pass

def readPairwise(indir, outfile):
	with open(outfile, "w") as output:
		# Open output file and write header
		output.write("TranscriptID,dN,dS,dN/dS,lnL\n")
		infiles = glob(indir + "*")
		for infile in infiles:
			# Get gene id and number of species remaining for each file
			filename = infile.split("/")[-1]
			transcript = filename.split(".")[0]
			with open(infile, "r") as mlc:
				for line in mlc:
					if line[:2] == "t=":
						# Extract dN, dS, and dN/dS
						i = line.index("dN/dS")
						j = line.index("dN ")
						k = line.index("dS ")
						dnds = line[i:j].split("=")[1].strip()
						dn = line[j:k].split("=")[1].strip()
						ds = line[k:].split("=")[1].strip()
					elif "lnL" in line:
						lnl = line.split("=")[1].strip()
			# Append data to list and convert into string
			data = [dn, ds, dnds, lnl]
			transcript += ","
			for i in data:
				transcript += str(i) + ","
			transcript = transcript[:-1]
			output.write(transcript + "\n")			

def main():
	parser = argparse.ArgumentParser(description="This script will \
concatenate CodeML output and print them in a single file.")
	parser.add_argument("--multiple", action="store_true",
 help="Indicates that a multiple alignment was used")
	parser.add_argument("--pairwise", action="store_true",
 help="Indicates that a pairwise alignment was used.")
	parser.add_argument("-i", help="Path to CodeML output directory")
	parser.add_argument("-o", help="Name of output file.")
	# Parse command
	args = parser.parse_args()
	indir = args.i
	if indir[-1] != "/":
		indir += "/"
	outfile = args.o
	multiple = args.multiple
	pairwise = args.pairwise
	concatenate(multiple, pairwise, indir, outfile)

if __name__ == "__main__":
	main()
