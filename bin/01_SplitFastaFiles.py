'''This will take the aligned multiple FASTA file with multiple genes 
    and create individual FASTA alignment files for each gene.


    Copyright 2016 by Shawn Rupp'''

from sys import argv

def splitFasta(infile, path):
    # Open input file and split into one alignment per gene
    with open(infile, "r") as fasta:
        newid = True
        seq = ""
        n = 0
        for line in fasta:
            if line != "\n":
                # Concatenate lines for each gene
                seq += str(line)
                if line[0] == ">":
                    # Determine number of sequences and species names
                    n += 1
                    if newid == True:
                        filename = str(line.split(".")[1])
                        newid = False
            elif line == "\n" and newid == False:
                if n >= 2:
                    # Print gene sequences to file if there are at least two
                    # species and reset for next gene
                    outfile = (path + "01_splitFastaFiles/" + filename + "."
                                + str(n) + ".fa")
                    with open(outfile, "w") as output:
                            output.write(seq)
                    newid = True
                    seq = ""
                    n = 0
                else:
                    pass

def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python 01_splitFastaFiles.py \
<input fasta alignment> <path to output directory>")
        quit()
    else:
        infile = argv[1]
        path = argv[2]
        splitFasta(infile, path)

if __name__ == "__main__":
    main()
