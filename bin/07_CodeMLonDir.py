'''This program will run CodeML on a directory of single gene alignments.
It will generate a unique control file and tree file for each input gene
before invoking CodeML using the number of CPUs specified by the user
(default =1).

    Copyright 2016 by Shawn Rupp'''

from sys import argv
from glob import glob
from subprocess import Popen
from shlex import split
from functools import partial
from multiprocessing import cpu_count
from multiprocessing.dummy import Pool
import re
import os

# Define max number of threads
MAXCPU = cpu_count()

def controlFiles(path):
    '''Reads input files and stores them in memory'''
    usertree = False
    with open(path + "codeml.ctl", "r") as control:
        ctl = control.readlines()
        for line in ctl:
            # Determine if a phylogenic tree is needed
            if "runmode = 0" in line or "runmode = 1" in line:
                usertree = True
    # Only run if a tree is required
    if usertree == True:
        with open(path + "codeml.tree", "r") as intree:
            tree = intree.readlines()[0].rstrip()
        # Create list of species and coresponding nodes
        nodes = {}
        splt = tree.split(",")
        for i in splt:
            i = i.replace("(", "")
            i = i.replace(")", "")
            i = i .replace(";", "")
            if "$" in i:
                print()
                print("\tPlease only use the pound sign (#) to indicate nodes.")
                print()
                quit()
            if "#" in i:
                node = i.split()
                nodes[node[0]] = node[1]
        with open(path + "ref.tree", "w") as reftree:
            # Add semicolon to the end of the tree if it is not present
            if tree[-1] != ";":
                tree += ";"
            # Remove node markers and write corrected tree to file
            if "#" in tree:
                for j in range(0, tree.count("#")):
                    i = tree.index("#")
                    tree = tree[:i-1] + tree[i+2:]
            reftree.write(tree)
    elif usertree == False:
        pass
    return usertree, nodes, ctl

def rmSeqs(path, retainStops):
    '''Creates a dictionary of sequences removed from specific genes'''
    rmseqs = {}
    # Only include sequences with internal stop codons if the have been removed
    if retainStops == False:
        with open(path + "internalStops.txt", "r") as stops:
            # Add sequences with internal stop codons
            for line in stops:
                splt = line.split("\t")
                if splt[0] == "Gene":
                    pass
                elif splt[0] in rmseqs:
                    rmseqs[splt[0]].append(str(splt[1]).rstrip())
                else:
                    # Must typecast dictionary values in order to create a list of
                    # strings
                    rmseqs[splt[0]] = [str(splt[1]).rstrip(), ]
    with open(path + "lowCount.txt", "r") as countlog:
        # Add sequences which were removed due to low content to the dictionary
        for line in countlog:
            splt = line.split("\t")
            if splt[0] == "Gene":
                pass
            if splt[0] in rmseqs:
                rmseqs[splt[0]].append(str(splt[1]).rstrip())
            else:
                rmseqs[splt[0]] = [str(splt[1]).rstrip(), ]
    return rmseqs

def runCodeml(usertree, path, ctl, rmseqs, nodes, file):
    '''Creates temporary control and tree files and runs CodeML.'''
    if file == "None":
        pass
    else:
        filename = file.split("/")[-1]
        geneid = filename.split(".")[0]
        # Set unique file names
        outfile = path + "07_codeml/" + geneid + ".mlc"
        tempctl = path + "07_codeml/" + geneid + ".ctl"
        tmpTree = path + "07_codeml/" + geneid + ".tree"
        if usertree == True:
            # Prune tree only if it will be used
            pruneTree(path, geneid, tmpTree, rmseqs, nodes)
        with open(tempctl, "w") as temp:
            # creates unique control file
            for line in ctl:
                if "seqfile" in line:
                    temp.write("\tseqfile = " + file + "\n")
                elif "outfile" in line:
                    temp.write("\toutfile = " + outfile + "\n")
                elif "treefile" in line:
                    if usertree == True:
                        temp.write("\ttreefile = " + tmpTree + "\n")
                else:
                    temp.write(line)
        # Calls CodeML
        cm = Popen(split("./paml/bin/codeml " + tempctl))
        cm.wait()
        os.remove(tempctl)
        os.remove(tmpTree)

def pruneTree(path, geneid, tmpTree, rmseqs, nodes):
    '''Calls the ape package in R to remove species whose sequences have been
removed due to low ccontent.'''
    # Identify removed sequences for file
    if geneid in rmseqs:
        # Call ape package to remove species from tree
        first = True
        outtree = path + "07_codeml/" + geneid + ".pruned"
        for i in rmseqs[geneid]:
            # Call R script
            if first == True:
                pt = Popen(split("Rscript bin/07_pruneTree.R " +
                            path + "ref.tree " + outtree + " " + i))
                pt.wait()
                first = False
            elif first == False:
                pt = Popen(split("Rscript bin/07_pruneTree.R " +
                            outtree + " " + outtree + " "  + i))
                pt.wait()
        # Read in ape output
        with open(outtree, "r") as pruned:
            tree = pruned.readlines()
        os.remove(outtree)
        # Add nodes back into tree
        for species in nodes:
            if species in tree:
                # Determine location and length of species name
                i = tree.index(species) + len(species)
                # Insert space and  node symbol after species name
                tree = (tree[:i] + " " + str(nodes[species]) +
                        tree[i:])
        with open(tmpTree, "w") as temptree:
            # Write tree to file
            string = ""
            for i in tree:
                string += i
            temptree.write(string)            
    else:
        # Read in ref tree instead of pruned tree
        with open(path + "ref.tree", "r") as ref:
            tree = ref.readlines()[0]
        for species in nodes:
            if species in tree:
                # Determine location and length of species name
                i = tree.index(species) + len(species)
                # Insert space and  node symbol after species name
                tree = (tree[:i] + " " + str(nodes[species]) +
                        tree[i:])
        with open(tmpTree, "w") as temptree:
            # Write tree to file
            string = ""
            for i in tree:
                string += i
            temptree.write(string)

def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python 07_CodeMLonDir.py <path to CodeML control file> \
-n <# of CPUs> <--retainStops>")
        quit()
    else:
        cpu = 1
        retainStops = False
        # Parse command
        for i in argv:
            if i == "-i":
                path = argv[argv.index(i) + 1]
            elif i == "-n":
                cpu = int(argv[argv.index(i) + 1])
            elif i == "--retainStops":
                retainStops = True
        # Set directory names and add a trailing "/" if necessary
        if path[-1] != "/":
            path += "/"
        # Make sure too many threads have not been specified
        if cpu > MAXCPU:
            cpu = MAXCPU
        # Reads in required data
        usertree, nodes, ctl = controlFiles(path)
        rmseqs = rmSeqs(path, retainStops)
        # Call CodeML for all files in a directory.
        inpath = path + "06_phylipFiles/" + "*.phylip"   
        files = glob(inpath)
        pool = Pool(processes = cpu)
        func = partial(runCodeml, usertree, path, ctl, rmseqs, nodes)
        pool.map(func, files)
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()
