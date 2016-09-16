'''This program will run KaKs_Calculator on a directory. Must be in the 
AlignmentProcessor directory to run.


    Copyright 2016 by Shawn Rupp'''

from sys import argv
from subprocess import Popen
from shlex import split
from glob import glob
import os

DEVNULL = open(os.devnull, "w")

def calculateKaKs(path):
    '''Calculates substition rates.'''
    inpath = path + "06_axtFiles/" + "*.axt"
    files = glob(inpath)
    for file in files:
        with open(file, "r") as infile:
            filename = file.split("/")[-1]
            # Create output file
            outfile = (path + "07_KaKsOutput/" + filename.split(".")[0]
                       + ".kaks")
            ck = Popen(split("bin/KaKs_Calculator -i " + file + " -o " +
                             outfile + " -m NG"), stdout = DEVNULL)
            ck.wait()

def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python 07_KaKsonDirectory.py \
<path to inut and output directories> <name of refernce species>")
        quit()
    else:
        path = argv[1]
        calculateKaKs(path)

if __name__ == "__main__":
    main()
