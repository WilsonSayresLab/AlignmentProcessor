'''This script will perform a permutation test on the means and medians of
a given column of one or two files.'''

import argparse
from datetime import datetime
import permute

def buildList1(header,column, value, data, infile):
	# Builds separate lists for condition 1 and condition 2
	with open(infile, "r") as indata:
		list1 = []
		list2 = []
		for line in indata:
			# Skip header
			if header == True:
				header = False
				pass
			elif header == False:
				# Evaluate conditional
				splt = line.split(",")
				cond = splt[column]
				# Append data to appropriate list and skip NAs
				if cond == value:
					try:
						list1.append(float(splt[data]))
					except ValueError:
						pass
				elif cond != value:
					try:	
						list2.append(float(splt[data]))
					except ValueError:
						pass
	return list1, list2

def buildList2(header, field1, field2, infile1, infile2):
	# Builds separate lists fr X and autosomal Ka, Ks, and Ka,KS values
	with open(infile1, "r") as file1:
		# Make list and tmeporary boolean
		list1 = []
		h = header
		for line in file1:
			# Skip header
			if h == True:
				h = False
				pass
			elif h == False:
				# Append values to list and skip NAs
				splt = line.split(",")
				try:
					list1.append(float(splt[field1]))
				except ValueError:
					pass
	with open(infile2, "r") as file2:
		# Make list and tmeporary boolean
		list2 = []
		h = header
		for line in file2:
			# Skip header
			if h == True:
				h = False
				pass
			elif h == False:
				# Append values to list and skip NAs
				splt = line.split(",")
				try:
					list2.append(float(splt[field2]))
				except ValueError:
					pass
	return list1, list2

def printPermuation(list1, list2, onesided, starttime):	   
	# Determine length of one list. The length of the other can be
	# determined by the difference between the given length and total length.
	length1 = len(list1)
	if onesided == False:
		print("\nPerforming two-sided permutation...")
		meancount, mediancount = permute.twoSided(list1, list2, length1)
		# Divides counts by number of replicates (10,000) to print p-values.
		print("\nP value for mean: ", meancount/10000)
		print("P value for median: ", mediancount/10000)
	elif onesided == True:
		print("\nPerforming one-sided permutation...")
		lmeancount, lmediancount, rmeancount, rmediancount = permute.oneSided(
														list1, list2, length1)
		# Divides counts by number of replicates (10,000) to print p-values.
		print("\nLeft-Tail (observed <= permuted average)") 
		print("P value for mean: ", lmeancount/10000)
		print("P value for median: ", lmediancount/10000)
		print("\nRight-Tail (observed >= permuted average)") 
		print("P value for mean: ", rmeancount/10000)
		print("P value for median: ", rmediancount/10000)
	# Prints elapsed time.
	print("\nTotal runtime: ", datetime.now() - starttime)	

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser(description = "This script will perform \
a permutation test on the means and medians of a given column of one or two \
files. Reuslts will be printed to the screen.")
	parser.add_argument("--noheader", action = "store_false", 
help = "Indicates that the file or files have no header.")
	parser.add_argument("--onesided", action = "store_true", 
help = "Indidcates that a one-sided permutation test should be conducted (compares both tails at once; default = two-sided).")
	parser.add_argument("--i1", help = "Path to first input file.")
	parser.add_argument("--i2", help = "Path to second input file (if using).")
	parser.add_argument("--c1", type = int,
help = "Column to permute from the first file.")
	parser.add_argument("--c2", type = int,
help = "column to permute from the second file (if using two files).")
	parser.add_argument("-c", help = "column from a single file that will be \
used to sort the data (e.g. 1=x).")
	args = parser.parse_args()
	# Determine number of input files to determine how to proceed
	if not args.i2:
		# Split sorting input
		c = args.c.split("=")
		column = c[0]
		value = c[1]
		list1, list2 = buildList1(args.noheader, column, value, args.c1, args.i1)
	elif args.i2:
		list1, list2 = buildList2(args.noheader, args.c1, args.c2, args.i1, args.i2)
	printPermuation(list1, list2, args.onesided, starttime)

if __name__ == "__main__":
	main()
