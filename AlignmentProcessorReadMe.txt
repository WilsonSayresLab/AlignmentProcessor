    Copyright 2016 by Shawn Rupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.


###############################################################
# AlignmentProcessor0.8 Package
#
#	Dependencies:	Python 3
#			Python 3 version of Biopython
#			Perl
#			PAML 
#			R
#			ape R package
###############################################################

### Contents ###
0. Introduction
1. Obtaining a fasta alignment
2. Running AlignmentProcessor
3. Individual Scripts
4. Outputs
5. Test

#-------------------------------
# 0. Introduction
#-------------------------------

AlignmentProcessor is a pipeline meant to quickly convert a multi-fasta 
alignment file into a format that can be read by KaKs_Calculator or PAML and 
optionally run those programs. When running codeml, alignment processor will
use input control and tree files as templates to create unique control and
tree files for each gene. It will call the ape R package to dynamically trim
the given phylogenic tree so that only species which remain in each gene's 
alignment after trimming are represented in that gene's tree. 

You can run the AlignmentProcessor wrapper which will call all of the python 
scripts in sequence, or you may call each script individually.

# Installing Biopython

Since Biopython offers the fastest and easiest ways to deal with alignments
in a python script, it is used in most of the scripts in this package.
Biopython is actually composed of pre-compiled Python scripts, so it is 
rather difficult to install directly. The easiest way is to download
a Python 3 version Anaconda (https://www.continuum.io/downloads). 
AlignmentProcessor was written in Python 3, but Python 3 will be installed 
alongside Anaconda, so you don't need to worry about that if you use 
Python 2.7. 

Once Anaconda is installed, all you have to do is paste the following into
a terminal and Anaconda will install Biopython for you:

	conda install -c https://conda.anaconda.org/anaconda biopython

# KaKs_Calculator

AlignmentProcessor0.8 is packaged with KaKs_Calculator2.0 binaries for Linux
and Windows, and a KaKs_Calculator1.2 binary for Mac (there is no 2.0 binary
available for OSX). Before using, copy or move the appropriate binary for your
system into the AlignmentProcessor bin which contains the python scipts.

# PAML 4.8

If you plan to use CodeML, you must first download PAML 
(http://abacus.gene.ucl.ac.uk/software/paml.html) and move the folder into the
AlignmentProcessor directory. Make sure that it is titled "paml".

# Ape

The most straightforward way to install ape, and most R packages, is through 
Bioconductor. If you do not have Bioconductor installed, open R and paste:

	source("https://bioconductor.org/biocLite.R")
	biocLite()

To install ape, enter:

	library("BiocInstaller")
	biocLite("ape")


#-------------------------------
# 1. Obtaining a fasta alignment
#-------------------------------
# UCSC Fasta Alignment
It is possible to download CDS fasta alignments from the UCSC Table browser.
This does, unfortunately, limit you to currently available alignments.
If they do have the alignment that you are interested in, however, it will 
probably be faster to download it rather than generate a new one. Since this
precludes the use of user-generated alignments, AlignmentProcessor has been 
written for Galaxy's Stitch Gene Blocks ouput. If you choose to use UCSC 
alignments, the sequence headers will have to be converted using the --ucsc 
option.

# User Generated Alignments
Since most alignments are in maf format, you will have to convert your 
alignment from maf to fasta. There seem to be very few programs that can do 
this; fortunately Galaxy's Stich Gene blocks not only converts a maf to fasta,
but it also separates sequences by genes, which is something we need to do 
anyway.

If you did not use a UCSC genome for the reference species in your alignment, 
you may need to upload the reference genome that you used as a custom build.
Make sure that the genome, maf file, and BED file are all set to the custom
build, and that the reference species build name is identical in all three 
files. Additionally, you may not be able to use the UCSC BED12 file if you did
not use a UCSC genome. If that is the case, you can either upload your own,
or, if you used an Ensembl genome, you can just remove the "chr_UN" and "chr"
chromosome prefixes from the file, and resubmit the file to Galaxy.

Upload your maf and BED files to Galaxy (usegalaxy.org) (or retrieve a BED 
file of the genes for your reference species using the UCSC Main link under 
the Get Data tab on the left of the screen). Select your species from the UCSC 
table browser, select Genes and Gene Predictions, ensembl genes, and the whole 
genome. Select BED as the output format and check the send to Galaxy box.

Once you have your data uploaded to Galaxy, select the Stich Gene Blocks tool 
under the Fetch Alignments/Sequences tab. Select your reference species' 
genome BED file in the Gene BED file dropdown menu. Change MAF Sourse to the 
maf file you uploaded (you may also use a locally catched alignment if it is 
available for all of your species of interest). Select the desired species 
IDs, leave Split into Gapless MAF Blocks to "no" (we will deal with gaps 
later on), and hit "Execute." When the tool has finished running, download the
output file.

This process will take a few hours, so plan accordingly.

Since this method offer the most flexibility for working with alignments, 
AlignmentProcessor was wrtien with this output format in mind and no further
formatting is required.

#-------------------------------
# 2. Running AlignmentProcessor
#-------------------------------

AlignmentProcessor is designed to convert the file into a useable format and 
run the substitutions quickly, so everything can be run with one command. Each
script can be run individually if necessary (each script's function and 
options will be discussed later). The input order for AlignmentProcessor's 
command line does not matter, but it does matter for the individual scripts.

To execute the AlignmentProcessor pipeline, you must first change into
the package directory. Otherwise it will not be able to locate the scripts
in the bin/ directory.

# Example Usage: 

	python AlignmentProcessor.py --ucsc --axt/phylip --kaks/codeml \
		--retainStops -% <decimal> -r <reference species> -i <input fasta file> \
		-o <path to output directory> 

# Required Arguments:

	-i	the path to your input fasta alignment.
	-o	the path to your working/output directory.
	-r	the build or common name of your reference species (more below).

# Optional Arguments:

	--ucsc	This will invoke 00_convertHeader.py, which will convert the 
		headers from UCSC fasta files so they only contain build
		names and gene IDs. This does not need to be run on Stich Gene
		Blocks output.

	--retainStops	This will tell AlignmentProcessor to retain sequences
			that contain internal stop codons. By default,
			sequences with internal stop codons will be removed 
			from the analysis as they may bias the results.

	--changeNames	Tells the program to change genome build names to 
			commom names (more below).

	-%	a decimal value specifying the minimum percentage of reads 
		that must remain after replacing unknown codons with gaps 
		(Default = 0.5). You may wish to use a lower threshold
		for highly diverged species or for low quality genomes.

	# Ka/Ks_Calculator/CodeML related arguments
	
	--axt/phylip	Specifies which format to convert the files into 
			(axt format for KaKs_Calculator; phylip for CodeML 
			or other programs; both commands may be run at once).

	--kaks	will run KaKs_Calculator (you much also specify --axt, 
		otherwise it will not run). Otherwise the program will quit 
		after converting the files.

	--codeml	will run codeml on all of the files in the 
			06_phylipFiles directory. You must also supply a 
			control file for CodeML which must be located in the 
			output directory you specified with the the -o option.
			This file must be titled "codeml.ctl" (the default 
			name given by PAML). 

	-n	if "--codeml" is selected, you may specify the number of CPUs 
		to run CodeML. CodeML itself cannot be parallelized, but 
		AlignmentProcessor can call multiple instances of CodeML to
		shorten overall run time. (Default = 1)

# Additional Commands
	
	-h/--help	will print the program's help dialogue

	-v/--version	will print the program version and copywright info

	--printNameList	will print the contents of 02_nameList.txt which 
			contains the list of genome builds and associated 
			common names. Must be run without any other arguments

	--addNameToList will add an entry to the 02_nameList.txt file.
			e.g. python AlignmentProcessor.py -- addNameToList \
				<build> <common name>
				

# Genome Builds and Common Names

	This specifies the reference species. To find it, either use the 
	--printNameList option or look in 02_nameList.txt in the bin. Check if
	the genome build is present in column 1 of the list. If it is, use the
	common name in column 2 of the list. If it is not, you may either use
	the --addNameToList option or add an entry to the file with the 
	build name that is present in your alignment as the first entry of a 
	new row, followed by a tab, then the desired common name. Make sure 
	there are no spaces in either name.

# The CodeML control file

	Codeml requires that all of its parameters be specified in one control 
	file (http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf). Provide a 
	control file with your desired parameters and AlignmentProcessor will 
	use it as template. It will only alter the input and output files so 
	that they are unique for each file. You may also need to provide a 
	tree file for codeml (see PAML manual above).

	The control file must be titled titled “codeml.ctl”, and it must be 
	located in the output directory.

# The CodeML tree file

	If your CodeML analysis requires a phylogneic tree, provide your 
	desired tree, titled "codeml.tree" in the output directory. Specify
	the tree as you want to appear to CodeML and be sure that the species
	names are specified as the common names in the 02_nameList.txt file, 
	but trimmed to ten characters (some programs still set a ten character
	limit on the length of names, so AlignmentProcessor trims the names). 
	The 07_CodeMLonDir.py script will save any nodes you have specified 
	with a "#" before sending a plain Newick tree to ape (which will not 
	work if there are PAML node symbols). It will then add any nodes back 
	into the tree after it has been trimmed. AlignmentProcessor will not 
	currently save nodes specified with "$" since it is difficult to 
	determine where a nested clade begins and ends.

# Invoking the Ka/Ks pipeline with a UCSC alignment:

	python AlignmentProcessor0.8.py --axt --kaks --ucsc -r green_anole \
	-i anolis_gallus.fa -o pairwiseKaKs/

# Invoking the CodeML pipeline with a de novo alignment:

	python AlignmentProcessor0.8.py --phylip --codeml -% 0.6 \
	-r green_anole -i anolis_gallus.fa -o codemlOutput/

#-------------------------------
# 3. Individual Scripts
#-------------------------------

Each script performs one or two functions on the input file or files, and 
saves the output to a new subdirectory. The location of the working directory
must be specified, but the subdirectories are hard-coded in the programs. 
Be sure to include the trailing "/" in the path. The AlignmentProcessor
wrapper will add it, but to avoid redundancies the individual scripts will 
not.

Remember that the order of the arguments does matter for these scripts.

# 00_ConvertHeader.py

	This script will convert headers for CDS fasta alignments from
	UCSC to be in the format: >"build_name"."gene_ID" 

	python convertHeader.py <path to input file>

# 01_SplitFastaFiles.py

	This script will split the input mult-fasta alignment into one file
	per gene. It will produce an output file for a gene if it has at least
	two sequences. 

	python 01_splitFastaFiles.py <input fasta alignment> \
		<path to output directory>

# 02_RemoveHeader.py

	This program will read through a directory that contains 
	aligned multiple FASTA files and replace FASTA headers with each 
	species' common name.

	python 02_RemoveHeaderOnDir.py 	<path to inut and output directories>

# 03_CheckFrame.py

	This script removes gaps introduced in the reference sequence by the
	alignment and removes corresponding sites in other species. It assumes
	that the reference sequnce was in frame before any gaps were inserted,
	and it returns the reference sequence to its original open reading
	frame. It will then replaces codons with missing nucleotides with gaps
	to remove unknown amino acids from the sequence.

	python 03_CheckFrameOnDir.py <path to inut and output directories> \
		<reference_species>

# 04_CountBases.py

	This program will check a multiple FASTA file to see that each species
 	retains a certain percentage of its nucleotide sequence. If not, it 
	will remove that sequence.

	Note: the AlignmentProcessor wrapper specifies 50% as a default value, 
	but the script itself does not, so you must specify one if you invoke 
	it on its own.

	python 05_CountBasesOnDir.py <threshold percentage as a decimal> \
		<path to inut and output directories>

# 05_ReplaceStopCodons.py

	This program will remove the internal stop codons (TAA, TAG, TGA)
	and replace with gaps (---) from the nucleotide alignment. Some 
	programs will not run properly if they enounter a premature stop
	codon.

	Terminal stop codons will be replaced, while sequences with internal 
	stop codons will have their gene id and sequence name recorded in the 
	internalStops.txt file. If --retainStops is specified, these sequences
	will be retained. If it is not, these sequences will be removed and 
	any which does not have at least two remaining sequences will not be 
	written to file.

	python 05_ReplaceStopCodonsOnDir.py \
		<path to inut and output directories> --retainStops(optional)

# 06_FASTAtoAXT.py

	This program executes 06_parseFastaIntoAXT.pl on an entire directory,
	allowing all of the contents of the directory to be converted to 
	axt files.

	Note: parseFastaIntoAXT.pl was provided by the developers of 
	KaKs_Calculator and, as such, is the only perl script in this package.

	python FASTAtoAXTonDirectory.py <path to inut and output directories>

# 06_FASTAtoPhylip.py

	This program will convert all files in an input directory
 	from fasta format to a phylip format.

	python 07_FASTAtoPhylip.py <number of species> \
		<path to inut and output directories>

# 07_KaKsonDir.py

	This program executes KaKs_Calculator on every file in a directory. 

	python 07_KaKsonDirectory.py <path to inut and output directories> \
		<name of refernce species>

# 07_CodeMLonDir.py

	This script will run codeml on every file in a directory. It requires
	the codeml.ctl file, and likely a tree file which it will supply to
	codeml. It will overwrite the "seqfile", "treefile", "outfile" lines 
	include the paths to the input phylip file, the output file, and the 
	tree file. It will also call the ape R package to trim the tree file 
	so that it only includes species which have not been filtered out.

	python 07_CodeMLonDir.py <path to codeml control file> \
		<path to input and output directories> \
		--retainStops(optional)

3 07_pruneTree.R

	This R script will call the ape package to dynamically trim input 
	trees for CodeML if any sequences have been removed. Species 
	whose sequences were removed in steps 4 or 5 will be removed from
	the temporary tree given to CodeML.

	(Caled by 07_CodeMLonDir.py)

# 08_compileKaKs_CSV.py

	This script concatonates the output from KaKs_Calculator into a text
	file. It adds a column for gene (or sequence) IDs, and prints the gene
	ID from the filename.

	python compileCSV.py <path to inut and output directories>

#-------------------------------
# 4. Outputs
#-------------------------------

A directory is created for each step, and, prior to the file conversion step, 
each directory will contain a series of single-gene fasta alignments. These 
are retained after the programs finish in case you need to examine them at any
point.

If --axt is selected, the 06_axtFiles directory will be populated with axt files 
for use with KaKs_Calculator. If --kaks is specified, KaKs_Calculator will run 
and its individual output files will be placed in the KaKsOutput directory. 
Finally, a single csv file containing all of output of KaKs_Calculator will be
printed to the output directory that you specified with the -o option.

If --phylip is specified, the program will convert the fasta files to phylip
after the stop codons have been removed. If you also specified --codeml, 
AlignmentProcessor will edit and submit the control file to CodeML and the 
output files will be saved in 07_codeml. Since there are many different things 
that can be done with the codeml output files, AlignmentProcessor does not 
attmept to concatenate specific parts from the output files.

If you wish to convert convert the files to both formats, specify both --axt
and --phylip the program will convert the fasta files to both formats. You may
also run one of the individual scripts on the 07_rmStops directory to convert
the files in a separate step. AlignmentProcessor will not, however, run 
KaKs_Calculator and CodeML simultaneously, as this could require too much 
memory.

#-------------------------------
# 5. Run AlignmentProcessor on test data
#-------------------------------

# To test KaKs_Calculator:
Change directory into the AlignmentProcessor folder. Paste the followig into
a terminal:

python AlignmentProcessor.py --axt --kaks --ucsc -r anoCar2 \
-i kaksTest.fa -o test/

This will return a text file with 11 lines.

# To test CodeML:
The test directory already contains sample CodeML control and tree files, so
all you need  to do is change into the AlignmentProcessor direcotry and paste
the following:

python AlignmentProcessor.py --phylip --codeml --ucsc -n 2 -r anoCar2 \
-i codemlTest.fa -o test/

There should be 8 .mlc files in the 07_codeml directory.
