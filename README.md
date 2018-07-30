    Copyright 2016 by Shawn Rupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

# AlignmentProcessor1.6 Package

# Dependencies:	
	Python 3
	Python 3 version of Biopython
	PAML (if using CodeML)
	PhyML (if using CodeML for multiple species alignments)


AlignmentProcessor is a pipeline meant to quickly convert a multi-fasta 
alignment file into a format that can be read by KaKs_Calculator or PAML and 
optionally run those programs. When running CodeML, alignment processor will
use an input control file as a template to create a unique control for each 
gene. It will call PhyML to create a unique phylogenic tree for each gene 
based off of its sequences. 

You can run the AlignmentProcessor wrapper which will call all of the python 
scripts in sequence, or you may call each script individually.

# Installation

### GitHub
AlignmentProcessor is freely available on GitHub and can either be downloaded directly from the site,
can cloned using the following command:

	git clone https://github.com/WilsonSayresLab/AlignmentProcessor.git

### pyfaidx  
AlignmentProcessor uses pyfaidx to parse fasta alignments. To install:

	pip install pyfaidx  

### KaKs_Calculator
AlignmentProcessor0.21 is packaged with KaKs_Calculator2.0 binaries for Linux
and Windows, and a KaKs_Calculator1.2 binary for Mac (there is no 2.0 binary
available for OSX). Before using, copy or move the appropriate binary for your
system into the AlignmentProcessor bin which contains the python scripts.

### PAML
If you plan to use CodeML, you must first download PAML 
(http://abacus.gene.ucl.ac.uk/software/paml.html) and move the folder into the
AlignmentProcessor directory. Make sure that it is titled "paml".

### PhyML
If you plan to use CodeML, you must also download PhyML 
(http://www.atgc-montpellier.fr/phyml/binaries.php). Similar to PAML, you must
move the folder into the AlignmentProcessor directory and change the name of
both the folder and the binary for your operating system to "PhyML".

# Quick Start

	python AlignmentProcessor.py --axt/phylip --kaks/codeml -r <reference species> -i <input fasta file> -o <path to output directory> 

### Please see AlignmentProcessorReadMe.pdf for specific options and instructions.
