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
    # awk -F'\t' '{split($1, arr, "_"); $1=arr[1]; print}' OFS='\t' braker_noiso.gtf > noiso.gtf # this might not work if the header ID has an underscore in it, for ex XP_19871203.1_Drosophila
    python ../../remove_isoforms/getAnnoFastaFromJoingenes.py -g ../*fa.masked -o "${BASE_NAME}"_proteome_noiso -f braker_noiso.gtf


    # Copy the braker.aa file to the main directory and rename it
    cp "${BASE_NAME}"_proteome_noiso.aa ../"${BASE_NAME}"_proteome_noiso.fa
    cd ../
    echo "we are currently in the $BASE_NAME directory:"
    pwd
    # Deactivating the previous conda environment and activating the correct one
    conda deactivate
    conda activate busco_annotation
    #################################################################################################

   
    busco -i "${BASE_NAME}_proteome_noiso.fa" -o ${BASE_NAME}_noiso_busco_annotated -m protein -l $LINEAGE --cpu $THREADS -f
    rm -r "busco_downloads"
    

    # Deactivating the previous conda environment to run interproscan in base
    conda deactivate
    #################################################################################################

    echo "Running Interproscan"
    sed -i 's/*//g' "${BASE_NAME}_proteome_noiso.fa"
    bash ~/my_interproscan/interproscan-5.69-101.0/interproscan.sh -i "${BASE_NAME}_proteome_noiso.fa" -f tsv -d "$BASE_NAME" -cpu $THREADS

    # Check if the braker.aa file was produced and is not empty
    if [ -f "${BASE_NAME}_proteome_noiso.fa.tsv" ] && [ -s "${BASE_NAME}_proteome_noiso.fa.tsv" ]; then
        echo "Interproscan ran successfully."
    else
        echo "Error: proteome.fa.tsv file not found or is empty in the directory $BASE_NAME, please recheck your run."
        exit 1
    fi

    #################################################################################################

    echo "Processing final proteome file"
    # Files
    input_file="${BASE_NAME}_proteome_noiso.fa.tsv"
    proteome_file="${BASE_NAME}_proteome_noiso.fa"
    header_file="${BASE_NAME}_annotated_headers_noiso.txt"
    output_file="${BASE_NAME}_annotated_proteome_noiso.fa"

    awk -F'\t' '
    {
        id = $1
        e_value = $9
        col6 = $6

        # Store the first line for each ID
        if (!(id in first_line)) {
            first_line[id] = ">" id " " col6
        }

        # Skip lines with "-" in e-value or column 6
        if (e_value == "-" || col6 == "-") {
            next
        }

        # Store the line with the lowest e-value for each ID
        if (!(id in min_e_value) || e_value < min_e_value[id]) {
            min_e_value[id] = e_value
            line[id] = ">" id " " col6
        }
    }
    END {
        for (id in first_line) {
            if (id in line) {
                print line[id]
            } else {
                print first_line[id]
            }
        }
    }' "$input_file" | sort -k1,1 > "$header_file"

    comm -13 \
        <(awk -F'\t' '{print ">" $1}' "$input_file" | sort -u) \
        <(grep ">" "$proteome_file" | sort -u) >> "$header_file"


    #######################################################

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
    }' $header_file > $temp_file

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

    # Clean up
    rm $temp_file

    # Check if the output file exists and is not empty
    if [ -s "$output_file" ]; then
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

}

remove_isoforms_meta(){
    BASE_NAME=$1
    cd $BASE_NAME
    remove_isoforms $BASE_NAME
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

        # Move all output files to the subdirectory
        mv *"noiso"* "$SUBDIR"
    done
}


source ~/miniconda3/etc/profile.d/conda.sh
THREADS=32
LINEAGE="euglenozoa"
main $1
