###############################################################
# This is a readme for analyzing Ka/Ks output from 
#	AlignmentProcessor
#
# 	Required programs:	R
#				Boot Package
###############################################################

#-------------------------------
# 1. Import Data into R
#-------------------------------

    trimmed<-read.csv("/media/yitzhak/Data1/Aves/KaKs1.2/galGal4KaKs.csv",header =TRUE)

	codeml <- read.csv("", header=TRUE)

#-------------------------------
# 2. Subset by chromosome
#-------------------------------

tz<-subset(trimmed,Chromosome.Name=="Z"&Chromosome.Name!="W",select=c(Ka,Ks,Ka.Ks))
tz<-na.omit(tz)
ta<-subset(trimmed,Chromosome.Name!="Z"&Chromosome.Name!="W",select=c(Ka,Ks,Ka.Ks))
ta<-na.omit(ta)

cz <- subset(codeml, Chromosome.Name=="Z"&Chromosome.Name!="W")
cz <- na.omit(cz)
ca <- subset(codeml, Chromosome.Name!="Z"&Chromosome.Name!="W")
ca <- na.omit(ca)

#-------------------------------
# 3. Generate histograms for Ka/Ks and Ks for the autosomes and Z chromosome
#-------------------------------

# Ka/Ks

{par(mfrow=c(2,1))
  hist(ta$Ka.Ks,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
         "G. gallus-F. peregrinus Autosomal Ka/Ks",xlab="Ka/Ks")
  hist(tz$Ka.Ks,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
         "G. gallus-F. peregrinus Z Chromosome Ka/Ks",xlab="Ka/Ks")}
# Ks

{par(mfrow=c(2,1))
 hist(ta$Ks,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
        "G. gallus-F. peregrinus Autosomal Ka/Ks",xlab="Ka/Ks")
 hist(tz$Ks,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
        "G. gallus-F. peregrinus Z Chromosome Ka/Ks",xlab="Ka/Ks")}

# dN/dS

{par(mfrow=c(2,1))
  hist(ca$dN.dS,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
         "G. gallus-F. peregrinus Autosomal dN/dS",xlab="dN/dS")
  hist(cz$dN.dS,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
         "G. gallus-F. peregrinus Z Chromosome dN/dS",xlab="dN/dS")}
# Ks

{par(mfrow=c(2,1))
 hist(ca$dS,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
        "G. gallus-F. peregrinus Autosomal dS",xlab="dS")
 hist(cz$dS,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
        "G. gallus-F. peregrinus Z Chromosome dS",xlab="dS")}

#-------------------------------
# 4. Determine mean Ka, Ks, and Ka/Ks
#-------------------------------
# Load bootstrap function and boot
source("/home/yitzhak/Dropbox/Scripts/R/bootstrap.R")

# Ka/Ks
# Autosomes
bootstrap(ta$Ka)
bootstrap(ta$Ks)
bootstrap(ta$Ka.Ks)
# Z chromosome
bootstrap(tz$Ka)
bootstrap(tz$Ks)
bootstrap(tz$Ka.Ks)

# dN/dS
# Autosomes
bootstrap(ca$dN)
bootstrap(ca$dS)
bootstrap(ca$dN.dS)
# Z chromosome
bootstrap(cz$dN)
bootstrap(cz$dS)
bootstrap(cz$dN.dS)

#-------------------------------
# 5. Compare mean Ka/Ks using Wilcoxon Tests between autosomes and Z
#-------------------------------

wilcox.test(tz$Ka.Ks,ta$Ka.Ks)
wilcox.test(cz$dN.dS,ca$dN.dS)
wilcox.test(trimmed$Ka.Ks,codeml$dN.dS)
