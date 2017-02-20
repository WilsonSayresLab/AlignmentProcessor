###############################################################
#-------------- OverView ----------------
#
#  This script will cal AlignmentProcessor with increasing percentage thresholds
#	on human to mouse, platypus, and aardvark pairwise alignments.
#
#	Required Programs:	Python3
#				Perl
#				Biopython
#				AlignmentProcessor0.9
#				KaKs_Calculator2.0
#############################################################

#-------------------------------
# 0. Download alignments and list of genes and loci
#-------------------------------

# Download alignments from UCSC

# Download a list of gene Ids, transcript IDs, and chromosome for hg19/GRCh37 from Ensembl:
# http://grch37.ensembl.org/biomart/martview/fb057aa8d2d5354e64fd94893ad5b13a

#-------------------------------
# 1. Make output directories
#-------------------------------

# Human to Mouse
	cd Mammal/Human_Mouse
	
	mkdir 0%
	mkdir 5%
	mkdir 10%
	mkdir 15%
	mkdir 20%
	mkdir 25%
	mkdir 30%
	mkdir 35%
	mkdir 40%
	mkdir 45%
	mkdir 50%
	mkdir 55%
	mkdir 60%
	mkdir 65%
	mkdir 70%
	mkdir 75%
	mkdir 80%
	mkdir 85%
	mkdir 90%
	mkdir 95%
	mkdir 100%

# Human to Aardvark

	cd Mammal/Human_Aardvark

	mkdir 0%
	mkdir 5%
	mkdir 10%
	mkdir 15%
	mkdir 20%
	mkdir 25%
	mkdir 30%
	mkdir 35%
	mkdir 40%
	mkdir 45%
	mkdir 50%
	mkdir 55%
	mkdir 60%
	mkdir 65%
	mkdir 70%
	mkdir 75%
	mkdir 80%
	mkdir 85%
	mkdir 90%
	mkdir 95%
	mkdir 100%

# Human to Platypus

	cd Mammal/Human_Platypus

	mkdir 0%
	mkdir 5%
	mkdir 10%
	mkdir 15%
	mkdir 20%
	mkdir 25%
	mkdir 30%
	mkdir 35%
	mkdir 40%
	mkdir 45%
	mkdir 50%
	mkdir 55%
	mkdir 60%
	mkdir 65%
	mkdir 70%
	mkdir 75%
	mkdir 80%
	mkdir 85%
	mkdir 90%
	mkdir 95%
	mkdir 100%

#-------------------------------
# 2. Change into AlignmentProcessor directory
#-------------------------------

	cd AlignmentProcessor/

#-------------------------------
# 2. Humann to Mouse
# Call AlignmentProcessor with for each percentage, sort output, join with orthologs file, and run R script for mean and medians
#-------------------------------

	echo "Percent,AutosomalMean,AutosomalMedian" > "Mammal/Human_Mouse/kaksAvg.tsv"

# 0%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.0 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/0%/

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 0%

# 5%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.05 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/5%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 5%

# 10%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.10 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/10%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 10%

# 15%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.15 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/15%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 15%

# 20%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.20 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/20%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 20%

# 25% 
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.25 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/25%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 25%

# 30%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.30 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/30%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 30%

# 35%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.35 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/35%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 35%

# 40%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.40 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/40%
	
	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 40%

# 45%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.45 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/45%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 45%

# 50%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.50 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/50%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 50%

# 55%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.55 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/55%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 55%

# 60%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.60 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/60%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 60%

# 65%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.65 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/65%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 65%

# 70% 
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.70 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/70%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 70%

# 75%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.75 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/75%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 75%

# 80%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.80 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/80%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 80%

# 85%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.85 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/85%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 85%

# 90%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.90 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/90%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 90%

# 95%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.95 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/95%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 95%

# 100%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 1.0 -i Mammal/Alignments/human_mouse.fa -o Mammal/Human_Mouse/100%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Mouse 100%

#-------------------------------
# 3. Humann to Aardvark
# Call AlignmentProcessor with for each percentage, sort output, join with orthologs file, and run R script for mean and medians
#-------------------------------

	echo "Percent,AutosomalMean,AutosomalMedian" > "Mammal/Human_Aardvark/kaksAvg.tsv"

# 0%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.0 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/0%/

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 0%

# 5%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.05 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/5%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 5%

# 10%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.10 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/10%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 10%

# 15%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.15 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/15%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 15%

# 20%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.20 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/20%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 20%

# 25% 
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.25 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/25%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 25%

# 30%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.30 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/30%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 30%

# 35%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.35 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/35%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 35%

# 40%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.40 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/40%
	
	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 40%

# 45%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.45 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/45%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 45%

# 50%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.50 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/50%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 50%

# 55%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.55 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/55%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 55%

# 60%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.60 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/60%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 60%

# 65%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.65 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/65%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 65%

# 70% 
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.70 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/70%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 70%

# 75%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.75 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/75%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 75%

# 80%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.80 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/80%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 80%

# 85%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.85 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/85%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 85%

# 90%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.90 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/90%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 90%

# 95%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.95 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/95%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 95%

# 100%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 1.0 -i Mammal/Alignments/Human_Aardvark.fa -o Mammal/Human_Aardvark/100%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Aardvark 100%

#-------------------------------
# 4. Humann to Platypus
# Call AlignmentProcessor with for each percentage, sort output, join with orthologs file, and run R script for mean and medians
#-------------------------------

	echo "Percent,AutosomalMean,AutosomalMedian" > "Mammal/Human_Platypus/kaksAvg.tsv"

# 0%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.0 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/0%/

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 0%

# 5%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.05 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/5%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 5%

# 10%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.10 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/10%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 10%

# 15%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.15 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/15%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 15%

# 20%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.20 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/20%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 20%

# 25% 
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.25 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/25%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 25%

# 30%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.30 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/30%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 30%

# 35%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.35 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/35%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 35%

# 40%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.40 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/40%
	
	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 40%

# 45%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.45 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/45%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 45%

# 50%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.50 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/50%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 50%

# 55%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.55 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/55%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 55%

# 60%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.60 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/60%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 60%

# 65%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.65 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/65%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 65%

# 70% 
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.70 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/70%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 70%

# 75%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.75 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/75%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 75%

# 80%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.80 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/80%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 80%

# 85%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.85 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/85%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 85%

# 90%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.90 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/90%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 90%

# 95%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 0.95 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/95%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 95%

# 100%
	python AlignmentProcessor.py --ucsc --axt --kaks -r Human -% 1.0 -i Mammal/Alignments/Human_Platypus.fa -o Mammal/Human_Platypus/100%

	Rscript ../ReadMe/kaksAverages.R Mammal/Human_Platypus 100%
