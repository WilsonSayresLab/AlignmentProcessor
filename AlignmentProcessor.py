'''AlignmentProcessor will run the subsituion rate pipeline to produce trimmed
axt or phylip files for use with KaKs_calculator or PhyMl.


    Copyright 2016 by Shawn Rupp

    This package is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.'''

from datetime import datetime
from sys import argv
from subprocess import Popen
from shlex import split
from glob import glob
import os

def makeDir(path, outdir, axt, phylip):
    '''Makes all sub-directories used by program.'''
    os.chdir(outdir)
    for i in ["01_splitFastaFiles", "02_rmHeader", "03_checkFrame",
              "04_countBasesPercent", "05_ReplaceStopCodons"]:
        try:
            os.mkdir(i)
        except FileExistsError:
            pass
    if axt == True:
        for i in ["06_axtFiles", "KaKsOutput"]:
            try:
                os.mkdir(i)
            except FileExistsError:
                pass
    if phylip == True:
        for i in ["06_phylipFiles", "07_codeml"]:
            try:
                os.mkdir(i)
            except FileExistsError:
                pass
    os.chdir(path)

def convertHeaders(fasta):
    '''Changes headers of UCSC CDS fasta files to includ eonly build name and
gene ID'''
    print("Converting fasta headers...")
    ch = Popen(split("python bin/00_ConvertHeader.py " + fasta))
    ch.wait()
    # Extcact converted file name
    path = fasta.split("/")[:-1]
    outpath = ""
    for i in path:
        outpath += i + "/"
    newfasta = outpath + "newHeader." + fasta.split("/")[-1]
    if ch.returncode == 0:
        return True, newfasta
            
def splitFasta(fasta, outdir):
    '''Splits fasta alignment into one file per gene.'''
    print("Splitting fasta file into separate files by gene...")
    sf = Popen(split("python bin/01_SplitFastaFiles.py " + " "
                + fasta + " " + outdir))
    sf.wait()
    if sf.returncode == 0:
        return True

def rmHeader(outdir):
    '''Removes header information and changes genome build to common name.'''
    print("Changing header...")
    rh = Popen(split("python bin/02_RemoveHeader.py " + outdir))
    rh.wait()
    if rh.returncode == 0:
        return True
    
def checkFrame(outdir, ref):
    '''Checks nucleotide frame of secondary species against the reference
species'''
    print("Checking nucleotide frames...")
    cf = Popen(split("python bin/03_CheckFrame.py " +  " " + outdir
                + " " + ref))
    cf.wait()
    if cf.returncode == 0:
        return True

def countBases(outdir, percent):
    '''Checks that each species has at least 50% of it's nucleotide content.'''
    print("Checking nucleotide content...")
    cb = Popen(split("python bin/04_CountBases.py " + " " + percent +
                     " " + outdir))
    cb.wait()
    if cb.returncode == 0:
        return True

def replaceStop(outdir):
    '''Removes stop codons for downstream analysis.'''
    print("Removing stop codons...")
    rs = Popen(split("python bin/05_ReplaceStopCodons.py " + " "
                     + outdir))
    rs.wait()
    if rs.returncode == 0:
        return True
    
def axtConvert(outdir, kaks, starttime):
    '''Open all input files in the directory and convert to axt script for
 use with KaKs_Calculator'''
    print("Converting fasta files into axt files...")
    ac = Popen(split("python bin/06_FASTAtoAXT.py " + outdir))
    ac.wait()
    if ac.returncode == 0:
        if kaks == True:
            return True
        elif kaks == False:
            print("Finished converting files.")
            elapsedtime = datetime.now() - starttime
            print("Total runtime: ", elapsedtime)

def calculateKaKs(outdir):
    '''Calculates substition rates.'''
    print("Calculating Ka/Ks values...")
    ck = Popen(split("python bin/07_KaKsonDir.py " + outdir))
    ck.wait()
    if ck.returncode == 0:
        return True

def printCSV(outdir, starttime):
    '''Prints Ka/Ks output as a single csv file.'''
    print("Printing Ka/Ks output as a single text file...")
    Popen(split("python bin/08_compileKaKs.py " + outdir))
    print("Finished calculating Ka/Ks values.")
    elapsedtime = datetime.now() - starttime
    print("Total runtime: ", elapsedtime)

def phylipConvert(outdir, starttime, codeml):
    '''Converts fasta files to phylip format'''
    print("Converting fasta files to phylip...")
    pc = Popen(split("python bin/06_FASTAtoPhylip.py " + " " + outdir))
    pc.wait()
    if pc.returncode == 0:
        if codeml == False:
            elapsedtime = datetime.now() - starttime
            print("Total runtime: ", elapsedtime)
        return True

def runcodeml(outdir, starttime):
    '''Runs codeml on a directory.'''
    print("Running codeml...")
    cm = Popen(split("python bin/07_CodeMLonDir.py " + outdir + "codeml.ctl "
                     + outdir))
    cm.wait()
    if cm.returncode == 0:
        elapsedtime = datetime.now() - starttime
        print("Total runtime: ", elapsedtime)

def helplist():
    print()
    print("### AlignmentProcessor will run the subsituion rate pipeline to \
produce trimmed axt files for use with KaKs_calculator. ###")
    print("    example usage: python AlignmentProcessor.py -% <decimal> \
--axt/phylip --kaks/codeml -i <input fasta file> -o \
<path to output directory> -r <reference species>")
    print("    --axt    Converts files to axt for use in KaKs_Calcuator.")
    print("    --phylip   Converts files to phylip for use in PhyML.")
    print("    --kaks     Runs KaKs_Calcuator if --axt is also specified")
    print("    --codeml     Runs codeml if --phylip is also specified")
    print("    --ucsc    converts headers of CDS fasta files obtained from \
the UCSC genome browser")
    print("    -r    the name of the reference species which will be used to \
determine the open reading frame")
    print("    -%    Sets the percentage cutoff for the countBases step (50% \
by default).")
    print("    -i    path to input file containing a multifasta alignment")
    print("    -o    AlignmentProcessor will use this as its working \
directory and print output to this directory")
    print("    -v    prints coyright and version information to the screen")
    print("### Please note that you must be in the substituionManager \
directory to run the program.###")
    print()

def main():
    starttime = datetime.now()
    # Set optional parameters to False:
    axt = False
    phylip = False
    kaks = False
    codeml = False
    conv = False
    percent = "0.5"
    # Extract values from command line:
    for i in argv:
        if i == "-h" or i == "--help":
            helplist()
            quit()
        elif i == "-v" or i == "--version":
            print("\nAlignmentProcessor Copyright 2016 by Shawn Rupp\n")
            print("This program comes with ABSOLUTELY NO WARRANTY")
            print("This is free software, and you are welcome to redistribute\
 it under certain conditions\n")
            quit()
        elif i == "-i":
            fasta = argv[argv.index(i) + 1]
        elif i == "-o":
            outdir = argv[argv.index(i) + 1]
            if outdir[-1] != "/":
                outdir = outdir + "/"
        elif i == "-r":
            ref = argv[argv.index(i) + 1]
        elif i == "-%":
            percent = argv[argv.index(i) + 1]
        elif i == "--axt":
            axt = True
        elif i == "--phylip":
            phylip = True
        elif i == "--kaks":
            kaks = True
        elif i == "--codeml":                
            codeml = True
        elif i == "--ucsc":
            conv = True
    # Save working directory to variable to call other scripts:
    path = os.getcwd()
    path = path + "/"
    # Run program if the user did not ask for help:
    if "-h" not in argv and "--help" not in argv:
        # Set checkpoint variables to False:
        ch = False
        sf = False
        rh= False
        cf = False
        cb = False
        rs = False
        ac = False
        ck = False
        pc = False
        # Begin pipeline
        makeDir(path, outdir, axt, phylip)
        if conv == True:
            ch, newfasta = convertHeaders(fasta)
            if ch == True:
                sf = splitFasta(newfasta, outdir)
        else:
            sf = splitFasta(fasta, outdir)
        if sf == True:
            rh = rmHeader(outdir)
        if rh == True:
            cf = checkFrame(outdir, ref)
        if cf == True:
            cb =countBases(outdir, percent)
        if cb == True:
            rs = replaceStop(outdir)
        # Optionally covert files to axt format:
        if axt == True:
            if rs == True:
                ac = axtConvert(outdir, kaks, starttime)
            # Run KaKs_Calculator:
            if kaks == True and codeml == False:
                if ac == True:
                    ck = calculateKaKs(outdir)
                if ck == True:
                    printCSV(outdir, starttime)
        # Optionally covert files to phylip format:
        if phylip == True:
            if rs == True:
                pc = phylipConvert(outdir, starttime, codeml)
            # Run codeml
            if codeml == True and kaks == False:
                if pc == True:
                    runcodeml(outdir, starttime)

if __name__ == "__main__":
    main()
