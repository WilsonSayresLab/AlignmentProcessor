'''This program will remove the internal stop codons (TAA, TAG, TGA)
	and replace with gaps (---) from the nucleotide alignment. 
	Typically, this would be done to calculate dN/dS.

Note: This assumes prior filtering has been done so that:
   -the reference sequence is a conserved open reading frame.
   -each nucleotide sequence is divisible by three evenly.


    Copyright 2016 by Shawn Rupp'''

from sys import argv
from glob import glob
from Bio import AlignIO
import os

def openFiles(path, retainStops):
    '''Open all input files in the directory'''
    internalStops = path + "internalStops.txt"
    with open(internalStops, "w") as stops:
        stops.write("Gene\tSpecies\n")
    inpath = path + "04_countBasesPercent/" + "*.countBases"
    files = glob(inpath)
    for file in files:
        with open(file, "r") as infile:
	    # Determine number of remaining sequences per file
            n = 0
            for line in infile:
                if line[0] == ">":
                    n += 1
        with open(file, "r") as infile:
            filename = file.split("/")[-1]
            geneid = filename.split(".")[0]
            # Create output file:
            outfile = (path + "05_ReplaceStopCodons/" + geneid
                       + "." + str(n) + ".rmStops")
            # Call functions
            seqs = seqDict(infile, n)
            count, newseqs = removeStops(n, geneid, seqs, internalStops,
                                     retainStops)
            writeSeqs(count, newseqs, outfile, retainStops)

def seqDict(infile, n):
    '''Converts fasta into separate sequence objects, determines sequence names,
    and creates dictionary entries for each set of codons. Also replaces
    terminal stop codons with gaps.'''
    seqs = {}
    alignment = AlignIO.parse(infile, "fasta", seq_count=n)
    try:
        for item in alignment:
            for i in range(0, n):
                codons = []
                seq = str(item[i].seq)
                species = str(item[i].id)
                for j in range(0, len(seq), 3):
                    codons.append(seq[j:j +3].upper())
                    j += 3
                # Replace terminal stop codons so the program can identify
                # remaining internal stops
                if (codons[-1] == "TAA" or codons[-1] == "TAG" or
                    codons[-1] == "TGA"):
                    codons[-1] = "---"
                seqs[species] = codons
        return seqs
    except ValueError:
        pass
     
def removeStops(n, geneid, seqs, internalStops, retainStops):
    '''This will replace internal stop codons with gaps, create a list
of genes with internal stop codons, and remove those sequences.'''
    count = n
    try:
        for species in seqs:
            record = True
            for idx,codon in enumerate(seqs[species]):
                if codon == "TAA" or codon == "TAG" or codon == "TGA":
                # Record genes with internal stop codons
                    if record == True:
                        with open(internalStops, "a") as stops:
                            stops.write(geneid + "\t" + species + "\n")
                        record = False
                        count -= 1
                    if retainStops == True:
                        seqs[species][idx] = "---"
        if retainStops == "False":
            if record == False:
                # Remove sequence from alignment
                seqs.pop(species)
        return count, seqs
    except TypeError:
        pass
    
def writeSeqs(count, newseqs, outfile, retainStops):
    '''Writes output files if there at least two sequences remaining.'''
    with open(outfile, "w") as output:
        if retainStops == "False":
            if count >= 2:
                # Write out sequence if at least 2 sequences remain
                for species in newseqs:
                    output.write(">" + str(species) + "\n")
                    seq = ""
                    for codon in newseqs[species]:
                        seq += str(codon)
                    output.write(seq + "\n")
        elif retainStops == "True":
            # Write out all sequences
            for species in newseqs:
                output.write(">" + str(species) + "\n")
                seq = ""
                for codon in newseqs[species]:
                    seq += str(codon)
                output.write(seq + "\n")

def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python 05_ReplaceStopCodons.py \
<path to inut and output directories> <True/False>")
        quit()
    else:
        path = argv[1]
        retainStops = argv[2]
        openFiles(path, retainStops)

if __name__ == "__main__":
    main()
