'''This program will run CodeML on a directory of single gene alignments.


    Copyright 2016 by Shawn Rupp'''

from sys import argv
from glob import glob
from subprocess import Popen
from shlex import split
import re
import os

def readIn(path):
    '''Reads input files and stores them in memory'''
    usertree = False
    with open(path + "codeml.ctl", "r") as control:
        ctl = control.readlines()
        for line in ctl:
            # Set tree file name to tmp tree
            if "treefile" in line:
                i = ctl.index(line)
                ctl.pop(i)
                newline = "\ttreefile = " + path + "07_codeml/tmp.tree \n"
                ctl.insert(i, newline)
            # Determine if a phylogenic tree is needed
            elif "runmode" in line:
                if "0" in line or "1" in line:
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
            if "#" in i or "$" in i:
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
            if "$" in tree:
                for j in range(0, tree.count("$")):
                    i = tree.index("$")
                    tree = tree[:i-1] + tree[i+2:]
            reftree.write(tree)
    elif usertree == False:
        pass
    with open(path + "internalStops.txt", "r") as stops:
        # Creates a dictionary of sequences removed from specific genes
        rmseqs = {}
        for line in stops:
            splt = line.split("\t")
            if splt[0] in rmseqs:
                rmseqs[splt[0]].append(str(splt[1]).rstrip())
            else:
                # Must typecast dictionary values in order to create a list of
                # strings
                rmseqs[splt[0]] = [str(splt[1]).rstrip(), ]
    return usertree, nodes, ctl, rmseqs

def runCodeML(usertree, path, ctl, rmseqs, nodes, retainStops):
    '''Runs CodeML on all files in a directory using temporary controlfiles.'''
    # Open all input files in the directory
    inpath = path + "06_phylipFiles/" + "*.phylip"   
    files = glob(inpath)   
    for file in files:
        filename = file.split("/")[-1]
        geneid = filename.split(".")[0]
        outfile = path + "07_codeml/" + geneid
        tempctl = path + "07_codeml/tmp.ctl"
        if retainStops == False and usertree == True:
            # Prune tree only if it will be used
            pruneTree(path, geneid, rmseqs, nodes)
        with open(tempctl, "w") as temp:
            # creates unique control file
            for line in ctl:
                if "seqfile" in line:
                    temp.write("\tseqfile = " + file + "\n")
                elif "outfile" in line:
                    temp.write("\toutfile = " + outfile + ".mlc\n")
                else:
                    temp.write(line)
        # Calls CodeML
        cm = Popen(split("./paml/bin/codeml " + tempctl))
        cm.wait()
    os.remove(tempctl)
    os.remove(path + "07_codeml/tmp.tree")

def pruneTree(path, geneid, rmseqs, nodes):
    '''Calls the ape package in R to remove species whose sequences have been
removed due to low ccontent.'''
    # Identify removed sequences for file
    if geneid in rmseqs:
        first = True
        for i in rmseqs[geneid]:
            # Call R script
            if first == True:
                pt = Popen(split("Rscript bin/07_pruneTree.R " + path + " " +
                            path + "ref.tree " + i))
                pt.wait()
                first = False
            elif first == False:
                pt = Popen(split("Rscript bin/07_pruneTree.R " + path + " " +
                            path + "pruned.tree " + i))
                pt.wait()
        with open(path + "pruned.tree", "r") as pruned:
            tree = pruned.readlines()
        for species in nodes:
            if species in tree:
                i = tree.index(species) + len(species)
                tree = (tree[:i] + " " + str(nodes[species]) +
                        tree[i:])
        with open(path + "07_codeml/tmp.tree", "w") as temptree:
            string = ""
            for i in tree:
                string += i
            temptree.write(string)            
    else:
        with open(path + "ref.tree", "r") as ref:
            tree = ref.readlines()[0]
        for species in nodes:
            if species in tree:
                i = tree.index(species) + len(species)
                tree = (tree[:i] + " " + str(nodes[species]) +
                        tree[i:])
        with open(path + "07_codeml/tmp.tree", "w") as temptree:
            string = ""
            for i in tree:
                string += i
            temptree.write(string)

def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python 07_CodeMLonDir.py <path to codeml control file>")
        quit()
    else:
        retainStops = False
        path = argv[1]
        # Set directory names and add a trailing "/" if necessary
        if path[-1] != "/":
            path += "/"
        try:
            if argv[2] == "--retainStops":
                retainStops = True
        except IndexError:
            pass
        usertree, nodes, ctl, rmseqs = readIn(path)
        runCodeML(usertree, path, ctl, rmseqs, nodes, retainStops)

if __name__ == "__main__":
    main()
