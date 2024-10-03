#!/bin/bash

# The purpose of this script is to remove the isoforms from proteomes. This script was written after the species_genome_annotation.sh script was done running, and isoforms were not accounted for. So this script gets the isoforms from the print_longest_isoform.py script,  then uses the already annotated proteome to generate the new proteome with no isoforms

species=$1

sudo chmod 777 braker3/
cd braker3
sudo ~/Desktop/LAU-ATCG-main/BRAKER-report/scripts/predictionAnalysis/selectSupportedSubsets.py --fullSupport FULLSUPPORT --anySupport ANYSUPPORT --noSupport NOSUPPORT  braker.gtf hintsfile.gff
python ~/Desktop/LAU-ATCG-main/print_longest_isoform.py ANYSUPPORT > braker_noiso.gtf
awk -F'\t' '{split($1, arr, "_Pro"); $1=arr[1]; print}' OFS='\t' braker_noiso.gtf > noiso.gtf
python /home/microbiology/miniconda3/envs/annotation/bin/getAnnoFastaFromJoingenes.py -g ../*_genomic.fna.masked -o $species -f noiso.gtf
grep ">" ../*annotated_proteome.fa > header_file.txt

# Create a temporary file to hold the replacement data
temp_file=$(mktemp)

# Prepare the header dictionary in awk format
awk '{
    if (substr($0, 1, 1) == ">") {
        key = substr($1, 2)
        $1 = ""
        desc[key] = substr($0, 2)
    }
} END {
    for (k in desc) {
        print k "\t" desc[k]
    }
}' header_file.txt > $temp_file
output_file=../"$species"_final.fa
proteome_file="$species".aa
# Replace headers in the fasta file with the descriptions
awk -v header_file="$temp_file" '
BEGIN {
    while ((getline < header_file) > 0) {
        header[$1] = $2
    }
}
{
    if (substr($0, 1, 1) == ">") {
        id = substr($1, 2)
        if (id in header) {
            print ">" id " " header[id]
        } else {
            print $0
        }
    } else {
        print $0
    }
}' $proteome_file > $output_file


