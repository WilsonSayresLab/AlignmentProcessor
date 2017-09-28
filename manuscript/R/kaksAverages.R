###############################################################
#-------------- OverView ----------------
#
#  This R script will subset Ka/Ks data and print means and medians
#   for autosomes and the X chromosome and append to a csv file.
#
#	Required Programs:	R
#
# usage: Rscript kaksAverages.R <path to input file> 
#
#############################################################

# Catch arguments passed to R
args <- commandArgs(trailingOnly=TRUE)
path <- args[1]
percent <- args[2]

# Compile output file name
outfile <- file.path(path, "kaksAvg.tsv")

# Read input file and define output file
infile <- read.table(file.path(path,percent,"KaKs.tsv"),header=TRUE)
kaks <- na.omit(infile$Ka.Ks)

# Write average values to output file

avg<-c(percent,mean(kaks),median(kaks))
output <- t(avg)
write.table(output, file=outfile, append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE)
