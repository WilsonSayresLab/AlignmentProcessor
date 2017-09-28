###############################################################
# This is a readme for analyzing Ka/Ks output from 
#	AlignmentProcessor at different filtering point
#
# 	Required programs:	R
#				Boot Package
###############################################################

#-------------------------------
# 1. Import Data into R
#-------------------------------

filtered<-read.csv("/media/yitzhak/Data1/archive/Aves/kaksFiltering/noStops/galGal4-filteredKaKs.csv",header =TRUE)
stops<-read.csv("/media/yitzhak/Data1/archive/Aves/kaksFiltering/withStops/galGal4-withStopsKaKs.csv",header =TRUE)
unfiltered<-read.csv("/media/yitzhak/Data1/archive/Aves/kaksFiltering/unfiltered/galGal4-unfilteredKaKs.csv",header =TRUE)

# Remove NAs
filtered <- na.omit(filtered$Ka.Ks)
stops <- na.omit(stops$Ka.Ks)
unfiltered <- na.omit(unfiltered$Ka.Ks)

#-------------------------------
# 2. Determine median Ka, Ks, and Ka/Ks
#-------------------------------
# Load bootstrap function and boot
source("/home/yitzhak/Dropbox/Scripts/R/bootstrap.R")

# Ka/Ks
bootstrap(filtered)
bootstrap(stops)
bootstrap(unfiltered)

#-------------------------------
# 3. Import Autosomal and Z chromosome data to get transcript totals
#-------------------------------

fa <- read.csv("/media/yitzhak/Data1/archive/Aves/kaksFiltering/noStops/filteredKaKs-Autosomes.csv", header=TRUE)
fa <- na.omit(fa$Ka.Ks)
fz <- read.csv("/media/yitzhak/Data1/archive/Aves/kaksFiltering/noStops/filteredKaKs-Z.csv", header=TRUE)
fz <- na.omit(fz$Ka.Ks)

sa <- read.csv("/media/yitzhak/Data1/archive/Aves/kaksFiltering/withStops/withStopsKaKs-Autosomes.csv", header=TRUE)
sa <- na.omit(sa$Ka.Ks)
sz <- read.csv("/media/yitzhak/Data1/archive/Aves/kaksFiltering/withStops/withStopsKaKs-Z.csv", header=TRUE)
sz <- na.omit(sz$Ka.Ks)

ua <- read.csv("/media/yitzhak/Data1/archive/Aves/kaksFiltering/unfiltered/unfilteredKaKs-Autosomes.csv", header=TRUE)
ua <- na.omit(ua$Ka.Ks)
uz <- read.csv("/media/yitzhak/Data1/archive/Aves/kaksFiltering/unfiltered/unfilteredKaKs-Z.csv", header=TRUE)
uz <- na.omit(uz$Ka.Ks)

#-------------------------------
# 4. Calculate Median Ka/Ks for Autosomal and Z chromosome
#-------------------------------

bootstrap(fa)
bootstrap(fz)
bootstrap(sa)
bootstrap(sz)
bootstrap(ua)
bootstrap(uz)
