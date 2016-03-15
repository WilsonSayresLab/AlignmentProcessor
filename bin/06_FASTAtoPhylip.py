'''This program will convert all files in an input directory
 from fasta format to a phylip format.


    AlignmentProcessor0.6 Copyright 2016 by Shawn Rupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License (GPL3.txt) for more details.'''

from Bio import AlignIO 
from sys import argv
from glob import glob

def phylipConvert(path):
    '''Convert all input files in the directory to phylip and write to
output directory'''
    inpath = path + "05_ReplaceStopCodons/" + "*.rmStops"
    files = glob(inpath)
    for file in files:
        with open(file, "r") as infile:
            filename = file.split("/")[-1]
            # Extract number of sequnces from file name
            n = int(filename.split(".")[1])
            # Create output file
            outfile = (path + "06_phylipFiles/" + filename.split(".")[0]
                       + "." + str(n) + ".phylip")
            # Parse fasta file
            alignment = AlignIO.parse(infile, "fasta", n)
            # Convert alignemnt and write to file
            AlignIO.write(alignment, outfile , "phylip-sequential")
        # Add additional space after species name so codeml can read the file
        with open(outfile, "r") as infile:
            seq = infile.readlines()
        with open(outfile, "w") as output:
            output.write(seq[0])
            for line in seq[1:n+2]:
                newline = line[0:9] + "  " + line[10:]
                output.write(newline)
                
def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python 06_FASTAtoPhylip.py \
<path to inut and output directories>")
        quit()
    else:
        path = argv[1]
        phylipConvert(path)

if __name__ == "__main__":
    main()
