#!/bin/bash

remove_isoforms(){
    BASE_NAME=$1

    sudo chmod 777 braker3/
    cd braker3
    echo "we are currently in the braker directory:"
    pwd
    conda activate remove_isoform

    sudo ../../remove_isoforms/selectSupportedSubsets.py --fullSupport FULLSUPPORT --anySupport ANYSUPPORT --noSupport NOSUPPORT braker.gtf hintsfile.gff
    python ../../remove_isoforms/print_longest_isoform.py ANYSUPPORT > braker_noiso.gtf
    awk -F'\t' '{split($1, arr, "_"); $1=arr[1]; print}' OFS='\t' braker_noiso.gtf > noiso.gtf # this might not work if the header ID has an underscore in it, for ex XP_19871203.1_Drosophila
    python ../../remove_isoforms/getAnnoFastaFromJoingenes.py -g ../*masked -o "${BASE_NAME}"_proteome_noiso -f noiso.gtf


    # Copy the braker.aa file to the main directory and rename it
    cp "${BASE_NAME}"_proteome_noiso.aa ../"${BASE_NAME}"_proteome_noiso.fa


    cd ../
    echo "we are currently in the $BASE_NAME directory:"
    pwd
    # Deactivating the previous conda environment and activating the correct one
    conda deactivate

    proteome_noiso="${BASE_NAME}_proteome_noiso.fa"
    tsv="${BASE_NAME}_proteome.fa.tsv"
    headers="${BASE_NAME}_annotated_headers.txt"
    
    filter_tsv_by_proteome ${proteome_noiso} ${tsv}
    update_proteome_headers ${proteome_noiso} ${headers}

   
    # Check if the output file exists and is not empty
    if [ -s "${proteome_noiso%.fa}_annotated.fa" ]; then
        echo "************************************************************"
        echo "*                                                          *"
        echo "*               Isoforms Removed Succesfully               *"
        echo "*                                                          *"
        echo "*                        Hooray!                           *"
        echo "*                                                          *"
        echo "************************************************************"
    else
        echo "Error: final annotated fasta file not created or is empty, please recheck your run"
        exit 1
    fi

    mv "${proteome_noiso%.fa}_annotated.fa" ../finalProteomes/${BASE_NAME}
    mv "${tsv%.tsv}_noiso.tsv" ../finalProteomes/${BASE_NAME}


}

#updates the isoform-free proteome headers with their annotational description
update_proteome_headers () {
    local proteome_file="$1"
    local headers_file="$2"
    local output_file="${proteome_file%.fa}_annotated.fa"

    # Create a temporary associative array to store header mappings
    declare -A headers_map

    # Read the headers file and populate the associative array
    while IFS= read -r line; do
        header_name=$(echo "$line" | awk '{print $1}' | sed 's/>//')
        headers_map["$header_name"]="$line"
    done < "$headers_file"

    # Process the proteome file
    while IFS= read -r line; do
        if [[ "$line" =~ ^\> ]]; then  # Corrected quoting
            # Extract the header name (without the '>' sign)
            header_name=$(echo "$line" | awk '{print $1}' | sed 's/>//')

            # Check if there's a replacement in headers_map
            if [[ -n "${headers_map[$header_name]}" ]]; then
                # Replace the line with the matched line from headers_map
                echo "${headers_map[$header_name]}"
            else
                # If no match, print the original line
                echo "$line"
            fi
        else
            # Print sequence lines as they are
            echo "$line"
        fi
    done < "$proteome_file" > "$output_file"

    echo "Updated proteome saved as $output_file"
}


#retains tsv entries that are found in isoform-free proteome
filter_tsv_by_proteome() {
    
    local tsv_file="$2"
    local proteome_file="$1"
    local output_file="${tsv_file%.tsv}_noiso.tsv"
    local temp_headers_file=$(mktemp)

    # Extract headers from the proteome file and save to a temporary file
    grep "^>" "$proteome_file" | cut -d' ' -f1 | cut -c2- > "$temp_headers_file"

    # Use awk to filter TSV lines where the first column matches a header
    awk -F'\t' 'NR==FNR { headers[$1]; next } $1 in headers' "$temp_headers_file" "$tsv_file" > "$output_file"

    # Clean up the temporary file
    rm "$temp_headers_file"

    echo "Filtered output written to $output_file"
}

main(){
    # Read each line of the file into the BASE_NAMES array
    mapfile -t BASE_NAMES < "$1"

    mkdir -p finalProteomes

    for BASE_NAME in "${BASE_NAMES[@]}"; do
        echo "Processing $BASE_NAME..."

        # Create a subdirectory for each base name in finalProteomes
        SUBDIR="finalProteomes/$BASE_NAME"
        mkdir -p "$SUBDIR"

        cd $BASE_NAME
        remove_isoforms $BASE_NAME
    done
}


source ~/miniforge3/etc/profile.d/conda.sh
main $1
