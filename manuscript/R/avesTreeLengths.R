###############################################################
# This is a readme for analyzing tree length and log-likelihood
# from CodeML
#
# 	Required programs:	R (3.3.1)
#				Boot Package
#       Bootstrap.R Script
###############################################################

#-------------------------------
# 1. Import Data into R  
#-------------------------------

# AlignmentProcessor0.12
	h0 <- read.csv("CodeML/avesOutputH0.csv", header=TRUE)
	h1 <- read.csv("CodeML/avesOutputH1.csv", header=TRUE)

# AlignmentProcessor1.1
	h0 <- read.csv("branchSpecific/avesNullOutput.csv", header=TRUE)
	h1 <- read.csv("branchSpecific/avesAltOutput.csv", header=TRUE)

#-------------------------------
# 2. Subset 
#-------------------------------
	
# Select genes with at least 7 sequences remaining
	
	h07 <- subset(h0, Species>=7, select=c(TreeLength,lnL))
	h17 <- subset(h1, Species>=7, select=c(TreeLength,lnL))

# Select genes with at least 14 sequences remaining

	h014 <- subset(h0, Species>=14, select=c(TreeLength,lnL))
	h114 <- subset(h1, Species>=14, select=c(TreeLength,lnL))

#-------------------------------
# 3. Determine mean and median tree lengths
#-------------------------------
# Load bootstrap function and boot
source("R/bootstrap.R")

# H0 Min. 3
bootstrap(h0$TreeLength)

# H0 Min. 7
bootstrap(h07$TreeLength)

# H0 Min. 14
bootstrap(h014$TreeLength)

# H1 Min. 3
bootstrap(h1$TreeLength)

# H1 Min. 7
bootstrap(h17$TreeLength)

# H1 Min. 14
bootstrap(h114$TreeLength)

#-------------------------------
# 4. Compare tree length means using Wilcoxon Tests
#-------------------------------

# Between models
wilcox.test(h0$TreeLength,h1$TreeLength)
wilcox.test(h07$TreeLength,h17$TreeLength)
wilcox.test(h014$TreeLength,h114$TreeLength)

# Within models
wilcox.test(h0$TreeLength,h07$TreeLength)
wilcox.test(h0$TreeLength,h014$TreeLength)
wilcox.test(h07$TreeLength,h014$TreeLength)

wilcox.test(h1$TreeLength,h17$TreeLength)
wilcox.test(h1$TreeLength,h114$TreeLength)
wilcox.test(h17$TreeLength,h114$TreeLength)

#-------------------------------
# 5. Plot log-likelihood 
#-------------------------------

{par(mfrow=c(2,1))
 hist(h0$lnL,breaks=seq(-250000,0,100),xlim=c(-20000,0),ylim=c(0,250),
      main="Aves Alignment Null Model Log-Likelihood",xlab="Log-Likelihood")
 hist(h1$lnL,breaks=seq(-230000,1000,100),xlim=c(-20000,0),ylim=c(0,250),
 	main="Aves Alignment Alternative Log-Likelihood",xlab="Log-Likelihood")
}

#-------------------------------
# 6. Calculate log-likelihood ratio (G)
#-------------------------------

G <- 2 * (h1$lnL - h0$lnL)

# Determine proportion of genes for which the alternative model is better
# Degrees of freedom = 1 , alpha = 0.05

# critical value found in Statistical Table A in The Analysis of Biological Data, p. 704
sig <- G >= 7.88
sum(sig)

# Calculate proportion of significant genes
sum(sig)/7917

#-------------------------------
# 7. Determine number of genes with a given number of species present
#-------------------------------

sum(h1$Species=="15")

sum(h1$Species=="14")

sum(h1$Species=="13")

sum(h1$Species=="12")

sum(h1$Species=="11")

sum(h1$Species=="10")

sum(h1$Species=="9")

sum(h1$Species=="8")

sum(h1$Species=="7")

sum(h1$Species=="6")

sum(h1$Species=="5")

sum(h1$Species=="4")

sum(h1$Species=="3")
