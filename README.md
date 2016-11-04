    Copyright 2016 by Shawn Rupp and Melissa Wilson Sayres

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.


##############################################################################
# AlignmentProcessor1.0 Package

#	Dependencies:	
			Python 3
			Python 3 version of Biopython
			PAML (if using CodeML)
			PhyML (if using CodeML)
##############################################################################


### Contents ###
0. Introduction
1. Obtaining a fasta alignment
2. Running AlignmentProcessor
3. Individual Scripts
4. Outputs
5. Test

-------------------------------
# 0. Introduction
-------------------------------

AlignmentProcessor is a pipeline meant to quickly convert a multi-fasta 
alignment file into a format that can be read by KaKs_Calculator or PAML and 
optionally run those programs. When running CodeML, alignment processor will
use an input control file as a template to create a unique control for each 
gene. It will call PhyML to create a unique phylogenic tree for each gene 
based off of its sequences. 

You can run the AlignmentProcessor wrapper which will call all of the python 
scripts in sequence, or you may call each script individually.

# Installing Biopython
Since Biopython offers the fastest and easiest ways to deal with alignments
in a python script, it is used in several of the scripts in this package.
The easiest way to install Biopython is to download a Python 3 version 
Anaconda (https://www.continuum.io/downloads). AlignmentProcessor was written
in Python 3, but Python 3 will be installed alongside Anaconda, so you don't 
need to worry about that if you use Python 2. 

Once Anaconda is installed, all you have to do is paste the following into
a terminal and Anaconda will install Biopython for you:

	conda install -c https://conda.anaconda.org/anaconda biopython

# KaKs_Calculator
AlignmentProcessor1.0 is packaged with KaKs_Calculator2.0 binaries for Linux
and Windows, and a KaKs_Calculator1.2 binary for Mac (there is no 2.0 binary
available for OSX). Before using, copy or move the appropriate binary for your
system into the AlignmentProcessor bin which contains the python scripts.

# PAML
If you plan to use CodeML, you must first download PAML 
(http://abacus.gene.ucl.ac.uk/software/paml.html) and move the folder into the
AlignmentProcessor directory. Make sure that it is titled "paml".

# PhyML
If you plan to use CodeML, you must also download PhyML 
(http://www.atgc-montpellier.fr/phyml/binaries.php). Similar to PAML, you must
move the folder into the AlignmentProcessor directory and change the name of
both the folder and the binary for your operating system to "PhyML".

-------------------------------
# 1. Obtaining a fasta alignment
-------------------------------

# UCSC Fasta Alignment
It is possible to download CDS fasta alignments from the UCSC Table browser.
This does, unfortunately, limit you to currently available alignments.
If they do have the alignment that you are interested in, however, it will 
be faster to download it rather than generate a new one. Since this
precludes the use of user-generated alignments, AlignmentProcessor has been 
written for Galaxy's Stitch Gene Blocks output. If you choose to use UCSC 
alignments, the sequence headers will have to be converted using the --ucsc 
option.

# User Generated Alignments
Since most alignments are in maf format, you will have to convert your 
alignment from maf to fasta. There seem to be very few programs that can do 
this; fortunately Galaxy's Stitch Gene blocks not only converts a maf to fasta,
but it also separates sequences by genes, which is something we need to do 
anyway.

If you did not use a UCSC genome for the reference species in your alignment, 
you may need to upload the reference genome that you used as a custom build.
Make sure that the genome, maf file, and BED file are all set to the custom
build, and that the reference species build name is identical in all three 
files. 

Additionally, you may not be able to use the UCSC BED12 file if you did
not use a UCSC genome. If that is the case, you can either upload your own,
or, if you used an Ensembl genome, you can just remove the "chr_UN" and "chr"
chromosome prefixes from the file, and resubmit the file to Galaxy. For NCBI 
genomes, you may download the gff from NCBI genome, use the UCSC utility 
"gff3ToGenePred" with the -useName and -honorStartStopCodons options, and 
use the UCSC utility "genePredToBed". This will return a BED12 file which 
may be submitted to Galaxy.

Upload your maf and BED files to Galaxy (usegalaxy.org) (or retrieve a BED 
file of the genes for your reference species using the UCSC Main link under 
the Get Data tab on the left of the screen). Select your species from the UCSC 
table browser, select Genes and Gene Predictions, ensembl genes, and the whole 
genome. Select BED as the output format and check the send to Galaxy box.

Once you have your data uploaded to Galaxy, select the Stitch Gene Blocks tool 
under the Fetch Alignments/Sequences tab. Select your reference species' 
genome BED file in the Gene BED file drop down menu. Change MAF Source to the 
maf file you uploaded (you may also use a locally catched alignment if it is 
available for all of your species of interest). Select the desired species 
IDs, leave Split into Gapless MAF Blocks to "no" (we will deal with gaps 
later on), and hit "Execute." When the tool has finished running, download the
output file.

This process will take a few hours, so plan accordingly.

Since this method offers the most flexibility for working with alignments, 
AlignmentProcessor was written with this output format in mind and no further
formatting is required.

-------------------------------
# 2. Running AlignmentProcessor
-------------------------------

AlignmentProcessor is designed to convert the file into a usable format and 
run the substitutions quickly, so everything can be run with one command. Each
script can be run individually if necessary (each script's function and 
options will be discussed later). The input order for AlignmentProcessor's 
command line does not matter, but it does matter for the individual scripts.

To execute the AlignmentProcessor pipeline, you must first change into
the package directory. Otherwise it will not be able to locate the scripts
in the bin/ directory.

If you are running CodeML and the program is interrupted, you may call the 
07_CodeMLonDir.py script and it will continue where CodeML left off. This
will save the time of having to run those genes through CodeML again (It will
do the same thing if you call the entire pipeline again, but there is no need 
to re-run the previous steps). It will not do the same for KaKs_Calculator 
since KaKs_Calculator is much faster, so it should not be a problem to just
invoke KaKs_Calculator on the whole directory again.

# Example Usage: 

	python AlignmentProcessor.py --ucsc --axt/phylip --kaks/codeml 
		--retainStops -% <decimal> -f <forward branch of codeml tree> 
		-r <reference species> -i <input fasta file> 
		-o <path to output directory> 

# Required Arguments:

	-i	the path to your input fasta alignment.
	-o	the path to your working/output directory.
	-r	the build or common name of your reference species (more below).

# Optional Arguments:

	--ucsc	This will invoke 00_convertHeader.py, which will convert the 
			headers from UCSC fasta files so they only contain build
			names and gene IDs. This does not need to be run on Stitch Gene
			Blocks output.

	--retainStops	This will tell AlignmentProcessor to retain sequences
					that contain internal stop codons. By default,
					sequences with internal stop codons will be removed 
					from the analysis as they may bias the results.

	--changeNames	Tells the program to change genome build names to 
					common names (more below).

	-%	a decimal value specifying the minimum percentage of reads 
		that must remain after replacing unknown codons with gaps 
		(Default = 0.5). You may wish to use a lower threshold
		for highly diverged species or for low quality genomes.

	
	--axt/phylip	Specifies which format to convert the files into 
					(axt format for KaKs_Calculator; phylip for CodeML 
					or other programs; both commands may be run at once).

	--kaks	will run KaKs_Calculator (you much also specify --axt, 
			otherwise it will not run). Otherwise the program will quit 
			after converting the files.

	-m	indicates the method for KaKs_Calculator to use to calculate 
		substitution rates (see below).

	--codeml	will run codeml on all of the files in the 
				06_phylipFiles directory. You must also supply a 
				control file for CodeML which must be located in the 
				output directory you specified with the the -o option.
				This file must be titled "codeml.ctl" (the default 
				name given by PAML). 

	-t	if "--codeml" is selected, you may specify the number of CPUs 
		to run CodeML. CodeML itself cannot be parallelized, but 
		AlignmentProcessor can call multiple instances of CodeML to
		shorten overall run time. (Default = 1)

	-f	the first ten characters (standard phylip format truncates the species
		names to ten characters) of the build name (or common name if you use 
		the --changeNames flag) of the species on the forward branch of the 
		phylogneic tree supplied by PhyML. This species does not have to be 
		the same as the reference species.

	--noCleanUp		tells the program to keep temporary CodeML control files,
					tree files, and other PhyML and CodeML output/temporary 
					files which will be located in the tmp directory. 
					These files are removed by default.

# Additional Commands
	
	-h/--help	will print the program's help dialogue

	-v/--version	will print the program version and copyright info

	--printNameList	will print the contents of 02_nameList.txt which 
			contains the list of genome builds and associated 
			common names. Must be run without any other arguments

	--addNameToList will add an entry to the 02_nameList.txt file.
			e.g. python AlignmentProcessor.py -- addNameToList 
				<build> <common name>
				

# Genome Builds and Common Names

	To find the common name of the refernece species, either use the 
	--printNameList option or look in 02_nameList.txt in the bin. Check if
	the genome build is present in column 1 of the list. If it is, use the
	common name in column 2 of the list. If it is not, you may either use
	the --addNameToList option or add an entry to the file with the 
	build name that is present in your alignment as the first entry of a 
	new row, followed by a tab, then the desired common name. Make sure 
	there are no spaces in either name.

# KaKs_Calculator Method

	KaKs_Calculator can calculate substitution rates in a number of different
	ways which can be specified to AlignmentProcessor using the "-m" flag. 
	These methods include estimations, that are generally much faster, and 
	maximum liklihood models, that should be more accurate. See the 
	KaKs_Calculator documentation in the KaKs_Calculotor folder for more 
	information.

	Estimations:
	NG (default in AlignmentProcessor)
	LWL
	LPB
	MLWL
	YN
	MYN
	
	Maximum Liklihood:
	GY
	MS (recommended fof maximum liklihood)

# The CodeML control file

	Codeml requires that all of its parameters be specified in one control 
	file (http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf). Provide a 
	control file with your desired parameters and AlignmentProcessor will 
	use it as template. 

	The control file must have a “.ctl” extension, and it must be 
	located in the output directory. Only provide one control file in this 
	directory. Examples are included with PAML and the controlFiles directory
	contains example control files for branch site, branch specific, and 
	pairwise analyses. These can simply be copied into your output directory 
	or you may supply one of your own.

# Invoking the Ka/Ks pipeline with a UCSC alignment:

	python AlignmentProcessor1.0.py --axt --kaks --ucsc -r anoCar2
	-i anolis_gallus.fa -o pairwiseKaKs/

# Invoking the CodeML pipeline with a de novo alignment:

	python AlignmentProcessor1.0.py --phylip --codeml -% 0.6
	-r anoCar2 -i anolis_gallus.fa -o codemlOutput/

-------------------------------
# 3. Individual Scripts
-------------------------------

Each script performs one or two functions on the input file or files and 
saves the output to a new subdirectory. The location of the working directory
must be specified, but the subdirectories are hard-coded in the programs. 
Be sure to include the trailing "/" in the path. The AlignmentProcessor
wrapper will add it, but to avoid redundancies the individual scripts will 
not.

Remember that the order of the arguments does matter for these scripts.

00_ConvertHeader.py

	This script will convert headers for CDS fasta alignments from
	UCSC to be in the format: >"build_name"."gene_ID" 

	python convertHeader.py <path to input file>

01_SplitFastaFiles.py

	This script will split the input multi-fasta alignment into one file
	per gene. It will produce an output file for a gene if it has at least
	two sequences. 

	python 01_splitFastaFiles.py <input fasta alignment>
		<path to output directory>

02_RemoveHeader.py

	This script will read through a directory that contains aligned multiple
	FASTA files and remove gene IDs from the fasta headers.	It will replace 
	build names with each species' common name if "--changeNames" is specified.

	python 02_RemoveHeaderOnDir.py <path to input and output directories>
		--changeNames(optional)

03_CheckFrame.py

	This script removes gaps introduced in the reference sequence by the
	alignment and removes corresponding sites in other species. It assumes
	that the reference sequence was in frame before any gaps were inserted,
	and it returns the reference sequence to its original open reading
	frame. It will then replace codons with missing nucleotides with gaps
	to remove unknown amino acids from the sequence.

	python 03_CheckFrameOnDir.py <path to input and output directories>
		<reference_species>

04_CountBases.py

	This program will check a multiple FASTA file to see that each species
 	retains a certain percentage of its nucleotide sequence. If not, it 
	will remove that sequence.

	Note: the AlignmentProcessor wrapper specifies 50% as a default value, 
	but the script itself does not, so you MUST specify one if you invoke 
	it on its own.

	python 05_CountBasesOnDir.py <threshold percentage as a decimal>
		<path to input and output directories>

05_ReplaceStopCodons.py

	This program will remove the internal stop codons (TAA, TAG, TGA)
	from the nucleotide alignment and replace them with gaps (---). Some 
	programs will not run properly if they encounter a premature stop codon.

	Terminal stop codons will be replaced, while sequences with internal 
	stop codons will have their gene id and sequence name recorded in the 
	internalStops.txt file. If --retainStops is specified, these sequences
	will be retained. If it is not, these sequences will be removed and 
	any gene which does not have at least two remaining sequences will not be
	written to file.

	python 05_ReplaceStopCodonsOnDir.py 
		<path to input and output directories> --retainStops(optional)

06_FASTAtoAXT.py

	This script converts the contents of a directory to axt format.

	python FASTAtoAXTonDirectory.py <path to input and output directories>

06_FASTAtoPhylip.py

	This program will convert all files in an input directory
 	from fasta format to a phylip format.

	python 07_FASTAtoPhylip.py <number of species> 
		<path to input and output directories>

07_KaKsonDir.py

	This program executes KaKs_Calculator on every file in a directory. 

	python 07_KaKsonDirectory.py <path to input and output directories>
		<method>

07_CodeMLonDir.py

	This script will run codeml on every file in a directory. It requires
	the codeml.ctl file. It will overwrite the "seqfile", "treefile", 
	"outfile" lines include the paths to the input phylip file, the output
	file, and the tree file. It will also run PhyML to create the tree file. If
	you are running CodeML and the program is interrupted, you may invoke this 
	script to pick up where you left off.
	
	Note: Since this script has greater utility as a stand-alone program, it 
	utilizes flags so that the order of the arguments does not matter. 

	python 07_CodeMLonDir.py -t <# of threads> -f <name of forward branch>
		-i <path to input and output directories> 


08_compileKaKs.py

	This script concatonates the output from KaKs_Calculator into a text
	file. It adds a column for gene (or sequence) IDs, and prints the gene
	ID from the filename.

	python compileKaKs.py <path to input and output directories>

-------------------------------
# 4. Outputs
-------------------------------

A directory is created for each step, and, prior to the file conversion step, 
each directory will contain a series of single-gene fasta alignments. These 
are retained after the programs finishes in case you need to examine them at
any point.

If --axt is selected, the 06_axtFiles directory will be populated with axt 
files for use with KaKs_Calculator. If --kaks is specified, KaKs_Calculator 
will run and its individual output files will be placed in the KaKsOutput 
directory. Finally, a single tsv file containing all of the output of 
KaKs_Calculator will be printed to the output directory that you specified 
with the -o option.

If --phylip is specified, the program will convert the fasta files to phylip
after the stop codons have been removed. If you also specified --codeml, 
AlignmentProcessor will edit and submit the control file to CodeML and the 
output files will be saved in 07_codeml. Since there are many different things 
that can be done with the codeml output files, AlignmentProcessor does not 
attempt to concatenate specific parts from the output files.

If you wish to convert the files to both formats, specifying both --axt
and --phylip will convert the fasta files to both formats. You may
also run one of the individual scripts on the 05_rmStops directory to convert
the files in a separate step. AlignmentProcessor will not, however, run 
KaKs_Calculator and CodeML simultaneously, as this could require too much 
memory.

-------------------------------
# 5. Run AlignmentProcessor on test data
-------------------------------

# To test KaKs_Calculator:
Change directory into the AlignmentProcessor folder. Paste the following into
a terminal:

python AlignmentProcessor.py --ucsc --axt --kaks -r anoCar2
-i test/kaksTest.fa -o test/

This will return a tsv file with 11 lines.

# To test CodeML:
The test directory already contains a sample CodeML control file, so
all you need  to do is change into the AlignmentProcessor directory and paste
the following:

python AlignmentProcessor.py --ucsc --phylip --codeml -t 2 -r anoCar2
-f anoCar2 -i test/codemlTest.fa -o test/

There should be 8 .mlc files in the 07_codeml directory.
