###############################################################
#-------------- OverView ----------------
#
#  This R script will call the ape package to dynamically trim input trees 
#	for CodeML if any sequences have been removed.
#
#	Required Programs:	R
#				ape
#
# usage: Rscript 07_pruneTree.R <path to input tree file> 
#   <path to output tree> <species to be excluded>
#
#	Copyright 2016 by Shawn Rupp
#############################################################

# Catch arguments passed to R
args <- commandArgs(trailingOnly=TRUE)

intree <- args[1]
outfile <- args[2]
exclude <- args[3]

# Split exclude into species
species <- strsplit(exclude, "-")[[1]]
l <- length(species)

# Load the ape package
library("ape")

# Read tree, prune extra branches, and write to file
tree <- read.tree(intree)
pruned <- drop.tip(tree, c(species[1:l]), trim.internal=FALSE)
write.tree(pruned, outfile)
