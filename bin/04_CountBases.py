'''This program will check a multiple FASTA file to see that each species
retains a certain percentage of its nucleotide sequence. If not, it will remove
that sequence.


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
from Bio import AlignIO
import os

def openFiles(percent, path):
    '''Opens all input files in the directory'''
    inpath = path + "03_checkFrame/" + "*.checkFrame"
    files = glob(inpath)
    for file in files:
        proceed = False
        with open(file, "r") as infile:
            filename = file.split("/")[-1]
            # Extract number of sequnces from file name
            n = int(filename.split(".")[1])
            # Create output file:
            outfile = (path + "04_countBasesPercent/" + filename.split(".")[0]
                       + "." + str(n) + ".countBases")
            proceed, seqs = seqDict(infile, n)
            if proceed == True:
                countBases(n, seqs, percent, outfile)

def seqDict(infile, n):
    '''Converts fasta into separate sequence objects, determine sequence names
    and create dictionary entries for each set of codons'''
    seqs = {}
    alignment = AlignIO.parse(infile, "fasta", seq_count=n)
    try:
        for item in alignment:
            for i in range(0, n):
                seq = str(item[i].seq)
                species = str(item[i].id)
                seqs[species] = seq
        return True, seqs
    except ValueError:
        return False, seqs

def countBases(n, seqs, percent, outfile):
    '''Counts the number of nucleotides and only writes the sequence to an
output file if they compose greater than the cutoff threshold of the sequence'''
    count = 0
    with open(outfile, "w") as output:
        for species in seqs:
            species = str(species)
            aligned = 0
            seq = seqs[species].upper()
            # Get count of all nucleotides
            gcount = seq.count("G")
            ccount = seq.count("C")
            acount = seq.count("A")
            tcount = seq.count("T")
            aligned += acount + tcount + ccount + gcount
            try:
                # Determine whether or not the sequence passes the treshold
                if aligned/len(seq) >= percent:
                    count += 1
                    output.write(">" + str(species) + "\n")
                    output.write(str(seqs[species]) + "\n")
            except ZeroDivisionError:
                pass
    # Delete output file if it does not have at least two sequences
    if count < 2:
        os.remove(outfile)
            

def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python 05_CountBasesOnDir.py \
<threshold percentage as a decimal> <path to inut and output directories>")
        quit()
    else:
        # Set directory names and add a trailing "/" if necessary:
        percent = eval(argv[1])
        path = argv[2]
        openFiles(percent, path)


if __name__ == "__main__":
    main()
