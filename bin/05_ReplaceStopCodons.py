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

def openFiles(path):
    '''Open all input files in the directory'''
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
            # Create output file:
            outfile = (path + "05_ReplaceStopCodons/" + filename.split(".")[0]
                       + "." + str(n) + ".rmStops")
            seqs = seqDict(infile, n)
            removeStops(seqs, outfile)

def seqDict(infile, n):
    '''Converts fasta into separate sequence objects, determine sequence names
    and create dictionary entries for each set of codons'''
    seqs = {}
    alignment = AlignIO.parse(infile, "fasta", seq_count=n)
    try:
        for item in alignment:
            for i in range(0, n):
                codons = []
                seq = str(item[i].seq)
                species = str(item[i].id)
                for j in range(0, len(seq), 3):
                    codons.append(seq[j:j +3])
                    j += 3
                seqs[species] = codons
        return seqs
    except ValueError:
        pass
     
def removeStops(seqs, outfile):
    '''This will remove stop codons from the sequences'''
    with open(outfile, "w") as output:
        try:
            for species in seqs:
                for idx,codon in enumerate(seqs[species]):
                    codon = codon.upper()
                    if codon == "TAA" or codon == "TAG" or codon == "TGA":
                        seqs[species][idx] = "---"
                output.write(">" + str(species) + "\n")
                seq = ""
                for codon in seqs[species]:
                    seq += str(codon)
                output.write(seq + "\n")
        except TypeError:
            pass
                        
def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python 05_ReplaceStopCodons.py \
<path to inut and output directories>")
        quit()
    else:
        path = argv[1]
        openFiles(path)

if __name__ == "__main__":
    main()
