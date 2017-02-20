#!/bin/bash

#SBATCH -n 14      
#SBATCH -N 1
#SBATCH -t 4-0:0
#SBATCH --job-name=AvesCodeML
#SBATCH -A mwilsons
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err  
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=smrupp@asu.edu

# Load modules:

module load python/3.5.1

# Change into program directory and copy control file to output folder

cd AlignmentProcessor/

cp controlFiles/branchSpecificNull.ctl /scratch/smrupp/aves/

#R un null model and concatenate results

python AlignmentProcessor.py --ucsc --phylip -r hg19 -t 14 -i /scratch/smrupp/aves/aves-hg19.fa -o /scratch/smrupp/aves/

python bin/ConcatenateCodeML.py --multiple -i /scratch/smrupp/aves/04_CodemlOutput/ -o /scratch/smrupp/aves/avesCodeMLNull.csv

# Rename output directory and replace control file

mv /scratch/smrupp/aves/04_CodemlOutput/*.mlc /scratch/smrupp/aves/avesH0/

mv /scratch/smrupp/aves/04_CodemlOutput/ /scratch/smrupp/aves/avesH1

rm /scratch/smrupp/aves/branchSpecificNull.ctl

rm /scratch/smrupp/aves/finishedCodeML.txt

cp controlFiles/branchSpecificAlt.ctl /scratch/smrupp/aves/

# Run alternative model

python bin/04_CallCodeML.py -t 14 -f falPer1 -i /scratch/smrupp/aves/03_phylipFiles/ -o /scratch/smrupp/aves/avesH1

python bin/ConcatenateCodeML.py --multiple -i /scratch/smrupp/aves/avesH1/ -o /scratch/smrupp/aves/avesCodeMLAlt.csv 
