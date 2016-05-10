'''This program will check a multiple FASTA file in a directory 
    to see that the open reading frame is conserved in each species,
    relative to reference species. If not it will correct it.
    This works only for insertions in non-reference species relative the
    reference. This program assumes that the reference species is in frame,
    so any non-triplet gaps in the reference will be removed from the
    alignments. This code is intended only for use with coding nucleotide
    alignments.

 For example it wil change this:
        hg19    ATT-TCATAG
        gorGor1 ATTTTCATAG
 To this:
        hg19    ATTTCATAG
        gorGor1 ATTTCATAG
 But this:
        hg 19   ATT---TCATAG
        gorGor1 ATTCTTTCATAG
 will not be removed because it doen't change the frame of the reference
 (hg19)


    Copyright 2016 by Shawn Rupp'''

from sys import argv
from glob import glob
from Bio import AlignIO
from collections import OrderedDict

def openFiles(path, ref):
    '''Opens all input files in the directory and creates output file names'''
    inpath = path + "02_rmHeader/" + "*.rmHeader"
    files = glob(inpath)
    for file in files:
        with open(file, "r") as infile:
            proceed = False
            filename = file.split("/")[-1]
            # Extract number of sequnces from file name
            n = int(filename.split(".")[1])
            # Create output file:
            outfile = (path + "03_checkFrame/" + filename.split(".")[0] + "." +
                       str(n) + ".checkFrame")
            proceed, seqs = seqDict(infile, n)
            if proceed == True:
                newseq = countFrame(seqs, ref)
                fixFrame(outfile, newseq)
            else:
                pass

def seqDict(infile, n):
    '''Convert fasta into separate sequence objects, determine sequence names
    and create dictionary entries for each set of codons'''
    seqs = OrderedDict()
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
        return True, seqs
    except ValueError:
        return False, seqs

def countFrame(seqs, ref):
    '''Removes gaps in the reference sequence introduced by
    alignment and removes corresponding sites in other species'''
    if ref in seqs:
        for codon in seqs[ref]: 
            if codon == "---":
                pass
            else:
                if "-" in str(codon):
                    i = codon.index("-")
                    idx = seqs[ref].index(codon)
                    codon = codon[:1] + codon[i + 1:]
                    for key in seqs:
                        triplet = seqs[key][idx]
                        triplet = triplet[:i] + triplet[i + 1:]
    return seqs
                        
def fixFrame(outfile, newseq):
    '''Replaces codons with missing nucleotides with gaps to 
remove unknown amino acids'''
    for species in  newseq:
        for codon in species:
            codon = codon.upper()
            if codon == "---":
                pass
            else:
                for i,char in enumerate(codon):
                    if char == "A" or char == "T" or char == "C" or char == "G":
                        pass
                    elif char == "-":
                        newseq[species][i] = "---"
                        break
    # Convert dictionary values to string and write to file
    with open(outfile, "w") as output:     
        for species in newseq:
            output.write(">" + str(species) + "\n")
            seq = ""
            for codon in newseq[species]:
                seq += str(codon)
            output.write(seq + "\n")

def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python 03_CheckFrame.py\
 <path to inut and output directories> <reference_species>")
        quit()
    else:
        path = argv[1]
        ref = argv[2]
        openFiles(path, ref)

if __name__ == "__main__":
    main()
