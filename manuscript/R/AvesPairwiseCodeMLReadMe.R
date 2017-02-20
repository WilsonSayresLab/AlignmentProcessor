###############################################################
# This is a readme for analyzing pairwise dN/dS output from 
#	CodeML
#
# 	Required programs:	R
#				Boot Package
###############################################################

#-------------------------------
# 1. Import Data into R and subset
#-------------------------------

input<-read.csv("galGal_falPer_dNdS.csv",header=TRUE,row.names=1)

# Subset by chromosome:
z<-subset(input,ChromosomeName=="Z",select=c(dN,dS,dNdS))
a<-subset(input,ChromosomeName!="Z"&ChromosomeName!="W",select=c(dN,dS,dNdS))

#-------------------------------
# 2. Generate histograms for dN/dS and dS for the autosomes and Z chromosome
#-------------------------------

# dN/dS

{par(mfrow=c(2,1))
  hist(a$dNdS,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
         "G. gallus-F. peregrinus Autosomal dN/dS",xlab="dN/dS")
  hist(z$dNdS,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
         "G. gallus-F. peregrinus Z Chromosome dN/dS",xlab="dN/dS")}
 
# dS

{par(mfrow=c(2,1))
 hist(a$dS,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
        "G. gallus-F. peregrinus Autosomal dN/dS",xlab="dN/dS")
 hist(z$dS,breaks=seq(0,21,0.01),xlim=c(0,1.5),ylim=c(0,750),main=
        "G. gallus-F. peregrinus Z Chromosome dN/dS",xlab="dN/dS")}

#-------------------------------
# 3. Determine mean dN, dS, and dN/dS
#-------------------------------

library("boot")

# Autosomes
mean(a$dN)
madN <- function(x,y) {return(mean(sample(a$dN, 300, replace=TRUE)))}
meanadN <- boot(data=a$dN,statistic=madN, R=1000)
ciadN <- boot.ci(meanadN, conf=0.95, type="norm")
ciadN$norm

mean(a$dS)
madS <- function(x,y) {return(mean(sample(a$dS, 300, replace=TRUE)))}
meanadS <- boot(data=a$dS,statistic=madS, R=1000)
ciadS <- boot.ci(meanadS, conf=0.95, type="norm")
ciadS$norm

mean(a$dNdS)
madNdS <- function(x,y) {return(mean(sample(a$dNdS, 300, replace=TRUE)))}
meanadNdS <- boot(data=a$dNdS,statistic=madNdS, R=1000)
ciadNdS <- boot.ci(meanadNdS, conf=0.95, type="norm")
ciadNdS$norm

# Z chromosome
mean(z$dN)
mzdN <- function(x,y) {return(mean(sample(z$dN, 300, replace=TRUE)))}
meanzdN <- boot(data=z$dN,statistic=mzdN, R=1000)
cizdN <- boot.ci(meanzdN, conf=0.95, type="norm")
cizdN$norm

mean(z$dS)
mzdS <- function(x,y) {return(mean(sample(z$dS, 300, replace=TRUE)))}
meanzdS <- boot(data=z$dS,statistic=madS, R=1000)
cizdS <- boot.ci(meanzdS, conf=0.95, type="norm")
cizdS$norm

mean(z$dNdS)
mzdNdS <- function(x,y) {return(mean(sample(z$dNdS, 300, replace=TRUE)))}
meanzdNdS <- boot(data=z$dNdS,statistic=mzdNdS, R=1000)
cizdNdS <- boot.ci(meanzdNdS, conf=0.95, type="norm")
cizdNdS$norm

#-------------------------------
# 4. Determine median dN, dS, and dN/dS
#-------------------------------

# Autosomes
median(a$dN)
madN <- function(x,y) {return(median(sample(a$dN, 300, replace=TRUE)))}
medianadN <- boot(data=a$dN,statistic=madN, R=1000)
ciadN <- boot.ci(medianadN, conf=0.95, type="norm")
ciadN$norm

median(a$dS)
madS <- function(x,y) {return(median(sample(a$dS, 300, replace=TRUE)))}
medianadS <- boot(data=a$dS,statistic=madS, R=1000)
ciadS <- boot.ci(medianadS, conf=0.95, type="norm")
ciadS$norm

median(a$dNdS)
madNdS <- function(x,y) {return(median(sample(a$dNdS, 300, replace=TRUE)))}
medianadNdS <- boot(data=a$dNdS,statistic=madNdS, R=1000)
ciadNdS <- boot.ci(medianadNdS, conf=0.95, type="norm")
ciadNdS$norm

# Z chromosome
median(z$dN)
mzdN <- function(x,y) {return(median(sample(z$dN, 300, replace=TRUE)))}
medianzdN <- boot(data=z$dN,statistic=mzdN, R=1000)
cizdN <- boot.ci(medianzdN, conf=0.95, type="norm")
cizdN$norm

median(z$dS)
mzdS <- function(x,y) {return(median(sample(z$dS, 300, replace=TRUE)))}
medianzdS <- boot(data=z$dS,statistic=madS, R=1000)
cizdS <- boot.ci(medianzdS, conf=0.95, type="norm")
cizdS$norm

median(z$dNdS)
mzdNdS <- function(x,y) {return(median(sample(z$dNdS, 300, replace=TRUE)))}
medianzdNdS <- boot(data=z$dNdS,statistic=mzdNdS, R=1000)
cizdNdS <- boot.ci(medianzdNdS, conf=0.95, type="norm")
cizdNdS$norm

#-------------------------------
# 5. Compare dN/dS distributions using Wilcoxon Tests
#-------------------------------

wilcox.test(a$dNdS,input$dNdS)
wilcox.test(z$dNdS,a$dNdS)

#-------------------------------
# 5. Compare dN/dS to Ks/Ks mean and medians with a permutation test
#-------------------------------

python permutation.py --c1 3 --c2 5 --i1 Chicken_Falcon_dNdS.csv --i2 galGal_falPerKaKs.trimmed.csv
