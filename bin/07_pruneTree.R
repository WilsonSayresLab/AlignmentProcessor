###############################################################
#-------------- OverView ----------------
#
#  This R script will call the ape package to dynamically trim input trees 
#	for CodeML if any sequences have been removed.
#
#	Required Programs:	R
#				ape
#
# usage: Rscript 07_pruneTree.R <path to output directory> 
#   <path to tree file> <species to be excluded>
#
#	Copyright 2016 by Shawn Rupp
#############################################################

# Catch arguments passed to R
args <- commandArgs(trailingOnly=TRUE)

path <- args[1]
tree <- args[2]
exclude <- args[3]

# Load the ape package
library(ape)

# Read tree, prune extra branches, and write to file
intree <- read.tree(tree)
prunedtree <- drop.tip(intree,exclude)
write.tree(prunedtree, file.path(path, "pruned.tree"))
