'''This script will concatenate CodeML output from AlignmentProcessor
 and print them in a single file.'''

import argparse
from glob import glob
import re

def inputList(indir):
	# Creates list of mlc files and determines which analysis was run
	infiles = glob(indir + "*")
	# Remove tmp dir
	for i in infiles:
		if i[i.rfind("/")+1:] == "tmp":
			infiles.remove(i)
	with open(infiles[0], "r") as ali:
		# Identify analysis
		for line in ali:
			if "dN/dS (w) for site classes" in line:
				typ = "brsite"
				break
			elif "dN & dS for each branch" in line:
				typ = "brspec"
				break
			elif "pairwise comparison, codon frequencies:" in line:
				typ = "pw"
				break
	return infiles, typ
		
def branchSite(infiles, outfile):
	# Saves data from branch site analysis to csv
	save = False
	with open(outfile, "w") as output:
		# Open output file and write header
		output.write("TranscriptID,Foreground_dN/dS,Background_dN/dS,\
TreeLength,lnL,Species,SitesUnderSelection (Amino Acid; PosteriorProbability)\n")
		for infile in infiles:
			try:
				# Get gene id and number of species remaining for each file
				filename = infile.split("/")[-1]
				transcript = filename.split(".")[0]
				species = filename.split(".")[1]
				if int(species) > 2:
					with open(infile, "r") as mlc:
						sites = ""
						for line in mlc:
							# Get substitution rates, tree lengths, etc
							if save == True:
								if "The grid" in line:
									save = False
								elif line.strip() and "Positive" not in line:
									# Record sites and probability
									splt = line.split()
									sites += ("{} ({}; {}),").format(splt[0],
															splt[1], splt[2])
							elif "tree length =" in line:
								length = line.split("=")[1].strip()
							elif "lnL" in line:
								lnl = line.split("):")[1]
								lnl = lnl.split()[0].strip()
							elif "background w" in line:
								# Convert multiple spaces to single space
								line = re.sub(" +", " ", line)
								bw = line.split()[4]
							elif "foreground w" in line:
								# Convert multiple spaces to single space
								line = re.sub(" +", " ", line)
								fw = line.split()[4]
							elif "Bayes Empirical Bayes" in line:
								save = True
					# Append data to list and convert into string
					data = [fw, bw, length, lnl, species, sites[:-1]]
					for i in data:
						transcript += "," + i
					output.write(transcript + "\n")
				else:
					# Skip files with only two sequences
					pass
			except IsADirectoryError:
				pass

def branchSpecific(infiles, outfile):
	# Saves data from branch specific analysis to csv
	with open(outfile, "w") as output:
		# Open output file and write header
		output.write("TranscriptID,dN,dS,dN/dS,TreeLength,lnL,Species\n")
		for infile in infiles:
			try:
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
					for i in data:
						transcript += "," + str(i)
					output.write(transcript + "\n")
				else:
					# Skip files with only two sequences
					pass
			except IsADirectoryError:
				pass

def pairwise(infiles, outfile):
	# Saves data from pairwise analysis to csv
	with open(outfile, "w") as output:
		# Open output file and write header
		output.write("TranscriptID,dN,dS,dN/dS,lnL\n")
		for infile in infiles:
			try:
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
			except IsADirectoryError:
				pass

def main():
	parser = argparse.ArgumentParser(description="This script will \
concatenate CodeML output and print them in a single file.")
	parser.add_argument("-i", help="Path to CodeML output directory")
	parser.add_argument("-o", help="Name of output file.")
	# Parse command
	args = parser.parse_args()
	indir = args.i
	if indir[-1] != "/":
		indir += "/"
	outfile = args.o
	infiles, typ = inputList(indir)
	if typ == "pw":
		pairwise(infiles, outfile)
	elif typ == "brsite":
		branchSite(infiles, outfile)
	elif typ == "brspec":
		branchSpecific(infiles, outfile)

if __name__ == "__main__":
	main()
