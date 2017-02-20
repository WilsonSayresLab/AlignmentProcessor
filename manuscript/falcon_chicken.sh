#!/bin/bash

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p smp
#SBATCH --mem-per-cpu=32310
#SBATCH -t 4-0:0
#SBATCH --job-name=gGal_fPer
#SBATCH -A mwilsons
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err  
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=smrupp@asu.edu

# Load the module for lastz
module load lastz/1.02.00

# Make working directory and copy genomes to it

mkdir lastz/galgal_falper/
cd lastz/galgal_falper/

cp /scratch/smrupp/genome/Gallus_gallus.Galgal4.dna_rm.toplevel.fa ./
cp /scratch/smrupp/genome/GCF_000337955.1_F_peregrinus_v1.0_genomic.fna ./

# Add path to lastz scripts

export PATH=$PATH:/home/smrupp/Whole_genome_alignment/pairwise/bin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/smrupp/Whole_genome_alignment/pairwise/bin/

# Invoke lastz pipeline:
lastz_CNM.pl Gallus_gallus.Galgal4.dna_rm.toplevel.fa GCF_000337955.1_F_peregrinus_v1.0_genomic.fna --run multi --cuts 100 --cpu 16 --hspthresh 2200 --inner 2000 --ydrop 3400 --gappedthresh 10000 --scores HoxD55 --chain --linearGap loose

# Change name and inline species identifiers of output maf file:
# mv all.maf galgal_falper.maf
# sed -i 's/target/galGal4./g' galgal_falper.maf
# sed -i 's/query/falPer1.0./g' galgal_falper.maf

