'''This program will read through a directory that contains 
    aligned multiple FASTA files and remove the NM identifier 
    (FASTA header) and all the trailing information (end lines)
    for each gene. It will also change all the 
    genome-build-names to common names.


    AlignmentProcessor0.6 Copyright 2016 by Shawn Rupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License (GPL3.txt) for more details.'''

from sys import argv
from glob import glob

def removeHeader(path):
    # Open all input files in the directory
    inpath = path + "01_splitFastaFiles/" + "*.fa"   
    files = glob(inpath)   
    for file in files:
        n = 0
        with open(file, "r") as infile:
            filename = file.split("/")[-1]
            # Create output file
            outfile = (path + "02_rmHeader/" + filename.split(".")[0] +
                       "." + filename.split(".")[1] + ".rmHeader")
            with open(outfile, "w") as output:
                for line in infile:
                    if line[0] == ">":
			# Split header to get build name and remove extra info
                        header = line.split(".")
                        name = changeNames(header[0])
			# Write common name to file
                        output.write(">" + str(name) + "\n")
                    else:
			# Write sequence to file in all caps
                        output.write(line.upper())

def changeNames(header):
    # Remove leading ">" from header to get the build name
    build = header[1:]
    speciesDict = {}
    with open("bin/02_nameList.txt", "r") as nameList:
        # Create dictionary of species names from file
        for line in nameList:
            speciesDict[line.split()[0]] = line.split()[1].rstrip()            
    # Compare the build name against the species dictionary. If it is matched
    # to a key, return the value (the common name).
        if build in speciesDict:
            name = speciesDict[build]
            return name

def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python 02_RemoveHeaderOnDir.py \
<path to inut and output directories>")
        quit()
    else:
        # Set directory names and add a trailing "/" if necessary
        path = argv[1]
        removeHeader(path)

if __name__ == "__main__":
    main()
