###############################################################
# This is a readme for analyzing pairwise fasta alignments
#	output from lastz.
#
# 	Required programs:	Lastz 
#				AlignmentProcessor1.2
#				PhyML
#				paml
###############################################################

#-------------------------------
# 1. Run Lastz for Alignment
#-------------------------------

# Download scripts from the Avian Genome Project and upload to server.
# Submit falcon_chicken.sh script on Saguaro. See comments at the end for changing file, target, and query names.

#-------------------------------
# 2. Run Stich Gene Blocks on Galaxy
#-------------------------------

# Upload the Lastz output maf file to Galaxy, and import the GalGal4 bed12 file from UCSC (select Genes and Gene Predictions,
# Ensembl genes, genome, and BED format before getting output. On the next page select whole gene). Remove the "chr_UN" and "chr"
# chromosome prefixes from the file, and resubmit the file to Galaxy. Also submit the chicken genome used in the alignment as a 
# custom build and set the builds of all three files to the custom build. (Alternatively, just use the UCSC chicken genome in the alignment.)

#-------------------------------
# 3. Run AlignmentProcessor
#-------------------------------

# KaKs_Calculator
	cd AlignmentProcessor/
	python AlignmentProcessor.py -t 4 --axt --kaks -r galGal4 -i Pairwise/galgal_falper.fa -o KaKs/

# CodeML 
	cd AlignmentProcessor/
	cp controlFiles/pairwiseAlt.ctl PairwiseCodeML1.2
	python AlignmentProcessor.py -t 4 --phylip --codeml -r galGal4 -i Pairwise/galgal_falper.fa -o PairwiseCodeML/

#-------------------------------
# 4. Prepare Output for analysis in R
#-------------------------------

# Download a list of gene and transcript IDs and their chromosome for galGal4 from Ensembl BioMart (v83). 
# Join this list with the Ka/Ks results to add gene, scaffold, and chromosome information to the results.
# Prior to joining, open each list in Excel and sort both files either ascending or descending by their transcript IDs.
# It does not matter which one you use, as long as you use the same one for each (join will not work properly if they are sorted differently).

# KaKs
	join -t "," --header --check-order -1 2 -2 1 "Pairwise/galGal4GeneTranscriptIDs.txt" "KaKs1.2/KaKs.csv" > "KaKs1.2/galGal4KaKs.csv"

# CodeML
	python bin/ConcatenateCodeML.py --pairwise -i PairwiseCodeML1.2/04_CodemlOutput -o PairwiseCodeML1.2/pairwiseAlt.csv
	join -t "," --header --check-order -1 2 -2 1 "Pairwise/galGal4GeneTranscriptIDs.txt" "PairwiseCodeML1.2/pairwiseAlt.csv" > "PairwiseCodeML1.2/galGal4Pairwise.csv"

#-------------------------------
# 5. Permutation test for Z chromosome vs. Autosomes
#-------------------------------

# Make two csv files from the joined file in step 4, one for only autosomal genes and one for only Z-linked genes. 
# Manually remove any genes with NAs for Ka, Ks, or Ka/Ks.

	cd PermutationScripts/
	
	python permutation.py --c1 7 --c2 7 -i1 KaKsOut/galGal_falPerAutosomalKaKs.csv -i2 KaKsOut/galGal_falPerZKaKs.csv


