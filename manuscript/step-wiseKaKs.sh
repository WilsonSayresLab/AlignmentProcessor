##############################################################################
# This script will call AlignmentProcessor on the chicken-falcon pairwise
# alignment and run KaKs_Calclator for each level of filtering.
#
#	Required:	AlignmentProcessor1.4
#				KaKs_Calculator
#
##############################################################################

#-----------------------------------------------------------------------------
# 1. Full filtering
#-----------------------------------------------------------------------------

python AlignmentProcessor.py -t 6 -r galGal4 --kaks -i /home/shawn/Documents/AvesAlignments/Pairwise/galgal_falper.fa -o /home/shawn/Documents/AvesAlignments/kaks/noStops

# Add ID and locus data
join -t "," --header --check-order -1 2 -2 1 "/home/shawn/Documents/AvesAlignments/Pairwise/galGal4GeneTranscriptIDs.csv" "/home/shawn/Documents/AvesAlignments/kaks/noStops/KaKs.csv" > "/home/shawn/Documents/AvesAlignments/kaks/noStops/galGal4-filteredKaKs.csv"

# Subset Autosomes and Z and run permutation script
cat "/home/shawn/Documents/AvesAlignments/kaks/noStops/galGal4-filteredKaKs.csv" | awk -F ',' '{if($3=="Z" || $3=="Chromosome Name") print $0}' > "/home/shawn/Documents/AvesAlignments/kaks/noStops/filteredKaKs-Z.csv"
cat "/home/shawn/Documents/AvesAlignments/kaks/noStops/galGal4-filteredKaKs.csv" | awk -F ',' '{if ($3 ~ /^[1-9]$/) print $0}' > "/home/shawn/Documents/AvesAlignments/kaks/noStops/filteredKaKs-Autosomes.csv"
cat "/home/shawn/Documents/AvesAlignments/kaks/noStops/galGal4-filteredKaKs.csv" | awk -F ',' '{if ($3 ~ /^[1-9][0-9]$/) print $0}' >> "/home/shawn/Documents/AvesAlignments/kaks/noStops/filteredKaKs-Autosomes.csv"

python permutation.py --c1 7 --c2 7 --i1 /home/shawn/Documents/AvesAlignments/kaks/noStops/filteredKaKs-Autosomes.csv --i2 /home/shawn/Documents/AvesAlignments/kaks/noStops/filteredKaKs-Z.csv

#-----------------------------------------------------------------------------
# 2. Retaining sequences with premature stops
#-----------------------------------------------------------------------------

python AlignmentProcessor.py -t 6 -r galGal4 --kaks --retainStops -i /home/shawn/Documents/AvesAlignments/Pairwise/galgal_falper.fa -o /home/shawn/Documents/AvesAlignments/kaks/withStops

# Add ID and locus data
join -t "," --header --check-order -1 2 -2 1 "/home/shawn/Documents/AvesAlignments/Pairwise/galGal4GeneTranscriptIDs.csv" "/home/shawn/Documents/AvesAlignments/kaks/withStops/KaKs.csv" > "/home/shawn/Documents/AvesAlignments/kaks/withStops/galGal4-withStopsKaKs.csv"

# Subset Autosomes and Z and run permutation script
cat "/home/shawn/Documents/AvesAlignments/kaks/withStops/galGal4-withStopsKaKs.csv" | awk -F ',' '{if($3=="Z" || $3=="Chromosome Name") print $0}' > "/home/shawn/Documents/AvesAlignments/kaks/withStops/withStopsKaKs-Z.csv"
cat "/home/shawn/Documents/AvesAlignments/kaks/withStops/galGal4-withStopsKaKs.csv" | awk -F ',' '{if ($3 ~ /^[1-9]$/ || $3=="Chromosome Name") print $0}' > "/home/shawn/Documents/AvesAlignments/kaks/withStops/withStopsKaKs-Autosomes.csv"
cat "/home/shawn/Documents/AvesAlignments/kaks/withStops/galGal4-withStopsKaKs.csv" | awk -F ',' '{if ($3 ~ /^[1-9][0-9]$/ || $3=="Chromosome Name") print $0}' >> "/home/shawn/Documents/AvesAlignments/kaks/withStops/withStopsKaKs-Autosomes.csv"

python permutation.py --c1 7 --c2 7 --i1 /home/shawn/Documents/AvesAlignments/kaks/withStops/withStopsKaKs-Autosomes.csv --i2 /home/shawn/Documents/AvesAlignments/kaks/withStops/withStopsKaKs-Z.csv

#-----------------------------------------------------------------------------
# 3. No filtering (using split fastas from previous steps)
#-----------------------------------------------------------------------------

python maskStops.py -i /home/shawn/Documents/AvesAlignments/kaks/noStops/01_splitFasta -o /home/shawn/Documents/AvesAlignments/kaks/unfiltered/rmStops

python bin/03_ConvertFasta.py --axt -i /home/shawn/Documents/AvesAlignments/kaks/unfiltered/rmStops -o /home/shawn/Documents/AvesAlignments/kaks/unfiltered/axt

python bin/04_CallKaKs.py -t 6 -i /home/shawn/Documents/AvesAlignments/kaks/unfiltered/axt -o /home/shawn/Documents/AvesAlignments/kaks/unfiltered/kaks

# Add ID and locus data
join -t "," --header --check-order -1 2 -2 1 "/home/shawn/Documents/AvesAlignments/Pairwise/galGal4GeneTranscriptIDs.csv" "/home/shawn/Documents/AvesAlignments/kaks/unfiltered/KaKs.csv" > "/home/shawn/Documents/AvesAlignments/kaks/unfiltered/galGal4-unfilteredKaKs.csv"

# Subset Autosomes and Z and run permutation script
cat "/home/shawn/Documents/AvesAlignments/kaks/unfiltered/galGal4-unfilteredKaKs.csv" | awk -F ',' '{if($3=="Z" || $3=="Chromosome Name") print $0}' > "/home/shawn/Documents/AvesAlignments/kaks/unfiltered/unfilteredKaKs-Z.csv"
cat "/home/shawn/Documents/AvesAlignments/kaks/unfiltered/galGal4-unfilteredKaKs.csv" | awk -F ',' '{if ($3 ~ /^[1-9]$/ || $3=="Chromosome Name") print $0}' > "/home/shawn/Documents/AvesAlignments/kaks/unfiltered/unfilteredKaKs-Autosomes.csv"
cat "/home/shawn/Documents/AvesAlignments/kaks/unfiltered/galGal4-unfilteredKaKs.csv" | awk -F ',' '{if ($3 ~ /^[1-9][0-9]$/ || $3=="Chromosome Name") print $0}' >> "/home/shawn/Documents/AvesAlignments/kaks/unfiltered/unfilteredKaKs-Autosomes.csv"

python permutation.py --c1 7 --c2 7 --i1 /home/shawn/Documents/AvesAlignments/kaks/unfiltered/unfilteredKaKs-Autosomes.csv --i2 /home/shawn/Documents/AvesAlignments/kaks/unfiltered/unfilteredKaKs-Z.csv

#-----------------------------------------------------------------------------
# 4. Permutations between Filtering levels
#-----------------------------------------------------------------------------

python permutation.py --c1 7 --c2 7 --i1 /home/shawn/Documents/AvesAlignments/kaks/unfiltered/galGal4-unfilteredKaKs.csv --i2 /home/shawn/Documents/AvesAlignments/kaks/withStops/galGal4-withStopsKaKs.csv

python permutation.py --c1 7 --c2 7 --i1 /home/shawn/Documents/AvesAlignments/kaks/unfiltered/galGal4-unfilteredKaKs.csv --i2 /home/shawn/Documents/AvesAlignments/kaks/noStops/galGal4-filteredKaKs.csv
