###############################################################
# This is a readme for analyzing Ka/Ks output from 
#	AlignmentProcessor.
#
# Output was returned using the KaKsMammals.sh script and kaksAverags.R.
#
# 	Required programs:	R
###############################################################

#-------------------------------
# 1. Import Data into R 
#-------------------------------

# Ka/Ks averages

mouse <- read.csv("Mammal/Human_Mouse/kaksAvg.tsv", header=TRUE)

aard <- read.csv("Mammal/Human_Aardvark/kaksAvg.tsv", header=TRUE)

plat <- read.csv("Mammal/Human_Platypus/kaksAvg.tsv", header=TRUE)

# Number of genes

genes <- read.csv("/home/yitzhak/Dropbox/AvesAlignments/Data/MammalGenesRemaining.csv", header=TRUE)

mgenes <- subset(genes, Species=="Mouse")
agenes <- subset(genes, Species=="Aardvark")
pgenes <- subset(genes, Species=="Platypus")

#-------------------------------
# 2. Create line plot for Ka/Ks values
#-------------------------------

{par(mfrow=c(3,2))
    
plot(mouse$Percent,mouse$AutosomalMean,type="o",col="green4",xlab="Threshold Percentage",ylab="dN/dS",
     pch=20,main="A. Human to Mouse Whole Genome\nMean and Median dN/dS",ylim=c(0,0.24))
    lines(mouse$Percent,mouse$AutosomalMedian,type="o",col="goldenrod3")
    legend("topright", c("Mean dN/dS","Median dN/dS"), col=c("green4","goldenrod3"), pch=c(20,21), lty=1, bty = "n")
plot(mgenes$Percent,mgenes$NumberOfTranscripts, type="o", xlab="Threshold Percentage", ylab="Number of Transcripts",
     main="B. Human to Mouse Trancripts Remaining after Filtering", col="blue3", ylim=c(0,13000))

plot(aard$Percent,aard$AutosomalMean,type="o",col="green4",xlab="Threshold Percentage",ylab="dN/dS",
     pch=20,main="C. Human to Aardvark Whole Genome\nMean and Median dN/dS",ylim=c(0,0.24))
    lines(aard$Percent,aard$AutosomalMedian,type="o",col="goldenrod3")
    legend("topright", c("Mean dN/dS","Median dN/dS"), col=c("green4","goldenrod3"), pch=c(20,21), lty=1, bty = "n")
plot(agenes$Percent,agenes$NumberOfTranscripts, type="o", xlab="Threshold Percentage", ylab="Number of Transcripts",
    main="D. Human to Aardvark Trancripts Remaining after Filtering", col="blue3", ylim=c(0,13000))

plot(plat$Percent,plat$AutosomalMean,type="o",col="green4",xlab="Threshold Percentage",ylab="dN/dS",
     pch=20,main="E. Human to Platypus Whole Genome\nMean and Median dN/dS",ylim=c(0,0.24))
    lines(plat$Percent,plat$AutosomalMedian,type="o",col="goldenrod3")
    legend("topright", c("Mean dN/dS","Median dN/dS"), col=c("green4","goldenrod3"), pch=c(20,21), lty=1, bty = "n")
plot(pgenes$Percent,pgenes$NumberOfTranscripts, type="o", xlab="Threshold Percentage", ylab="Number of Transcripts",
     main="F. Human to Platypus Trancripts Remaining after Filtering", col="blue3", ylim=c(0,13000))
}
