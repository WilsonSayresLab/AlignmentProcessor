###############################################################
# This is a readme for analyzing fasta alignmnets
#	downloaded from UCSC.
#
# 	Required programs:	AlignmentProcessor1.2
#						CodeML4.9
###############################################################

#-------------------------------
# 1. Alignments
#-------------------------------

# From the UCSC table browser page:
Select the human hg19 genome, genes and gene predictions, and ensembl genes.
Select CDS FASTA in ouput format, enter a file name, and select gzip compressed before hitting get output. On the next page,
Select the multiz100way alignment under MAF table. Select every species under birds, and deselect everything else. 

# From Ensembl BioMart:
# http://dec2013.archive.ensembl.org/biomart/martview/00c5ddbaf1b0544dc54db6cd39534529
Download a list of hg19/GRCh37 gene IDs with corresponding transcript IDs and chromosome. Save as a csv.

#-------------------------------
# 2. Run AlignmentProcessor to obtain pyhilp output
#-------------------------------

# Make sure to clone AlignmentProcessor from GitHub: https://github.com/WilsonSayresLab/AlignmentProcessor

Upload fasta alignment and avesCodeml.sh to cluster and submit job.

#-------------------------------
# 4. Add gene IDs and locus data 
#-------------------------------

# First download concatenated CodeML results from the cluster.
# Sort ascending by transcript ID (I did it in Excel since I inspected the files as well), then add gene IDs and locus data

	join -t "," --header --check-order -1 2 -2 1 "h19GeneTranscriptIDs.csv" "branchSpecific/avesCodeMLNull.csv" > "branchSpecific/avesNullOutput.csv"
	join -t "," --header --check-order -1 2 -2 1 "h19GeneTranscriptIDs.csv" "branchSpecific/avesCodeMLAlt.csv" > "branchSpecific/avesAltOutput.csv" 

#-------------------------------
# 5. Permute result files to obtain a p-value for mean and median tree lengths
#-------------------------------

# To compile the Cython scripts, change into the PermutationScripts directory and type: python setup.py build-ext --inplace

	cd PermutationScripts/

# dN
	python permutation.py --c1 3 --c2 3 --i1 branchSpecific/avesNullOutput.csv --i2 branchSpecific/avesAltOutput.csv

# dS
	python permutation.py --c1 4 --c2 4 --i1 branchSpecific/avesNullOutput.csv --i2 branchSpecific/avesAltOutput.csv

# dN/dS
	python permutation.py --c1 5 --c2 5 --i1 branchSpecific/avesNullOutput.csv --i2 branchSpecific/avesAltOutput.csv

# Tree length
	python permutation.py --c1 6 --c2 6 --i1 branchSpecific/avesNullOutput.csv --i2 branchSpecific/avesAltOutput.csv
