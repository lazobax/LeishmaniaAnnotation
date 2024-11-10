#!/bin/bash

# Function to process each genome accession
process_accession() {
    acc="$1"

    # Source conda environment
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate ncbi-datasets

    # Download dataset
    if ! datasets download genome accession "$acc" --include "protein"; then
        echo "$acc: Failed to download dataset" >> download_errors.txt
        return 1
    fi

    # Unzip the dataset
    if ! unzip "ncbi_dataset.zip"; then
        echo "$acc: Failed to unzip dataset" >> download_errors.txt
        # Move to zipFiles folder if unzipping fails
        mv "ncbi_dataset.zip" "zipFiles/${acc}.zip"
        return 1
    fi
    conda deactivate

    # Create proteomes directory if it doesn't exist
    mkdir -p proteomes

    # Move the protein file to the proteomes folder
    mv "ncbi_datasets/data/$acc/protein.faa" "proteomes/$acc.faa"

    # Activate BUSCO environment and run BUSCO
    conda activate busco
    busco -i "proteomes/$acc.faa" -o "busco_annotated_$acc" -m protein -l euglenozoa --cpu 32 -f

    # Extract the BUSCO summary and append to report.txt
    echo "$acc" >> "report.txt"
    grep "C:" "busco_annotated_$acc"/*txt >> "report.txt"
    
    conda deactivate
}

# List of genome accessions
accession_list=(
    "GCA_962240455.1"
    "GCA_000002845.2"
    "GCA_000227135.2"
    "GCA_000755165.1"
    "GCA_962240445.1"
    "GCA_962240435.1"
    "GCA_962239345.1"
    "GCA_000002875.2"
    "GCA_000234665.4"
    "GCA_017916335.1"
    "GCA_017916325.1"
    "GCA_017916305.1"
    "GCA_017918215.1"
    "GCA_017918225.1"
    "GCA_000002725.2"
)

# Loop through each accession and process it
for acc in "${accession_list[@]}"; do
    process_accession "$acc"
done
