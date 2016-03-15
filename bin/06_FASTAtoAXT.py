'''This program executes parseFastaIntoAXT.pl on an entire direcotry,
allowing all of the contents of the directory to be converted to axt files.


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
from subprocess import Popen
from shlex import split
from glob import glob
import os

def axtConvert(path):
    '''Open all input files in the directory and convert to axt script for
 use with KaKs_Calculator'''
    inpath = path + "05_ReplaceStopCodons/" + "*.rmStops"
    files = glob(inpath)
    for file in files:
        with open(file, "r") as infile:
            # Create output file:
            filename = file.split("/")[-1]
            outfile = (path + "06_axtFiles/" + filename + ".axt")
            ac = Popen(split("perl bin/06_parseFastaIntoAXT.pl " + file))
            ac.wait()
            # Moves file to output directory
            os.rename(file + ".axt",  outfile)

def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python 06_FASTAtoAXT.py \
<path to inut and output directories>")
        quit()
    else:
        path = argv[1]
        axtConvert(path)

if __name__ == "__main__":
    main()
