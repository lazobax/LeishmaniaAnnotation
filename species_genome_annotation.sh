#!/bin/bash

# Function to display usage message
usage() {
    echo "Usage: $0 -a <accession_number> -r <reference_fasta> -l <lineage> -o <output_directory> -t <threads>"
    echo "Options:"
    echo "  -a <accession_number>    NCBI Species accession number (must start with GCA_ or GCF_)"
    echo "  -r <reference_fasta>     Reference FASTA file (.fa, .fna, .fasta)"
    echo "  -l <lineage>             BUSCO analysis lineage term (https://busco.ezlab.org/list_of_lineages.html)"
    echo "  -o <output_directory>    Output directory (default: current working directory)"
    echo "  -t <threads>             Number of threads to use (default: 8)"
    echo "  -h                       Display this help message and exit"
    exit 1
}

# Initialize variables
ACCESSION=""
REFERENCE=""
LINEAGE=""
OUTPUT_DIR=$(pwd)  # Default output directory to the current working directory
THREADS=8  # Default value for threads

# Parse command line options using getopts
while getopts "a:r:l:o:t:h" opt; do
    case ${opt} in
        a )
            ACCESSION=$OPTARG
            ;;
        r )
            REFERENCE=$OPTARG
            ;;
        l )
            LINEAGE=$OPTARG
            ;;
        o )
            OUTPUT_DIR=$OPTARG
            ;;
        t )
            THREADS=$OPTARG
            # Check if THREADS is numeric
            if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
                echo "Error: Threads (-t) must be a numeric value."
                exit 1
            fi
            if [ $THREADS -gt $(nproc) ]; then
                echo "Error: You specified more threads than available ($THREADS > $(nproc)). Limited to the maximum available."
                THREADS=$(nproc)
            fi
            ;;
        h )
            usage
            ;;
    esac
done

# Check if required arguments are provided
if [ -z "$ACCESSION" ] || [ -z "$REFERENCE" ] || [ -z "$LINEAGE" ]; then
    usage
fi

# Check if the reference file exists
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference file $REFERENCE not found."
    exit 1
fi

# Check if the reference file has a valid extension (.fa, .fna, .fasta)
if ! [[ "$REFERENCE" =~ \.(fa|fna|fasta)$ ]]; then
    echo "Error: Reference file $REFERENCE does not have a valid .fa, .fna, or .fasta extension."
    exit 1
fi

# Determine if accession is GenBank or RefSeq
if [[ $ACCESSION =~ ^GCA_ ]]; then
    SOURCE="genbank"
elif [[ $ACCESSION =~ ^GCF_ ]]; then
    SOURCE="refseq"
else
    echo "Error: Accession number must start with GCA_ (GenBank) or GCF_ (RefSeq)"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run the ncbi-genome-download command with the appropriate source, flat-output option, and retry option
ncbi-genome-download invertebrate -s $SOURCE --assembly-accessions $ACCESSION --formats fasta --flat-output --output $OUTPUT_DIR -P -r 3

# Check if the download was successful
if [ $? -eq 0 ]; then
    echo "Download completed successfully and saved in $OUTPUT_DIR"
else
    echo "Error: Download failed"
    exit 1
fi

# Find the downloaded file
DOWNLOADED_FILE=$(find $OUTPUT_DIR -maxdepth 1 -type f -name "*.fna.gz")

# Extract the file if it exists
if [ -f "$DOWNLOADED_FILE" ]; then
    # Extract the file
    gunzip "$DOWNLOADED_FILE"

    # Get the base name of the file without the .fna.gz extension
    BASE_NAME=$(basename "$DOWNLOADED_FILE" .fna.gz)

    # Set the species FASTA file path
    SPECIES_FASTA="$OUTPUT_DIR/$BASE_NAME/$BASE_NAME.fna"

    # Create a directory with the base name
    mkdir -p "$OUTPUT_DIR/$BASE_NAME"
echo "$OUTPUT_DIR/$BASE_NAME"
chmod 777 "$OUTPUT_DIR/$BASE_NAME" 

    # Move the extracted .fna file to the new directory
    mv "$OUTPUT_DIR/$BASE_NAME.fna" "$SPECIES_FASTA"
else
    echo "Error: Downloaded file not found."
    exit 1
fi

# Change working directory to the newly created directory
cd "$OUTPUT_DIR/$BASE_NAME/"

# File to store the start and end times
log_file="$OUTPUT_DIR/$BASE_NAME/time_log.txt"

# Printing time and args to stdout and to logfile
start_time=$(date "+%Y-%m-%d %H:%M:%S")
start_seconds=$(date +%s)

echo "Start Time: $start_time"
echo "Start Time: $start_time" > "$log_file"
echo "Accession number: $ACCESSION"
echo "Accession number: $ACCESSION" >> "$log_file"
echo "Output directory: $OUTPUT_DIR/$BASE_NAME/"
echo "Output directory: $OUTPUT_DIR/$BASE_NAME/" >> "$log_file"
echo "Threads: $THREADS"
echo "Threads: $THREADS" >> "$log_file"

#################################################################################################

# Activating conda environment "annotation"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate annotation

# Building repeat database
BuildDatabase -name "$BASE_NAME" "$SPECIES_FASTA"

# Check for expected database files
DATABASE_FILES=("$BASE_NAME.translation" "$BASE_NAME.nhr" "$BASE_NAME.nin" "$BASE_NAME.njs" "$BASE_NAME.nog" "$BASE_NAME.nsq" "$BASE_NAME.nnd" "$BASE_NAME.nni")

for FILE in "${DATABASE_FILES[@]}"; do
    if [ ! -f "$FILE" ]; then
        echo "Error: Expected database file $FILE not found. Exiting."
        exit 1
    fi
done

# Run RepeatModeler with the specified number of threads
RepeatModeler -database "$BASE_NAME" -threads $THREADS

# Check if the families file was produced and is not empty
if [ -f "${BASE_NAME}-families.fa" ] && [ -s "${BASE_NAME}-families.fa" ]; then
    echo "Repeats successfully modeled."
else
    echo "Error: ${BASE_NAME}-families.fa was not produced or is empty. Exiting."
    exit 1
fi

# Run RepeatMasker using the produced families file
RepeatMasker -lib "${BASE_NAME}-families.fa" "$SPECIES_FASTA" -pa $THREADS
# -qq for faster but less accurate search, use for testing to speed up

# Check if the masked file was produced and is not empty
if [ -f "${SPECIES_FASTA}.masked" ] && [ -s "${SPECIES_FASTA}.masked" ]; then
    echo "Repeats masked successfully. The masked file is in the directory $OUTPUT_DIR"
else
    echo "Error: No .masked files found in the directory $OUTPUT_DIR, repeats not masked, please recheck your run."
    exit 1
fi

#################################################################################################

# Set Docker-related variables
file_name=$(basename "$SPECIES_FASTA")   # Set file_name as BASE_NAME
refdir=$(dirname "$REFERENCE")
ref_name=$(basename "$REFERENCE")
dirpath="$OUTPUT_DIR/$BASE_NAME"

# Echo the variables
echo "File name: $file_name"
echo "Reference directory: $refdir"
echo "Reference name: $ref_name"
echo "Output directory path: $dirpath"

# Run the Docker command for Braker3
docker run --user 1000:100 --rm -it \
 -v "${dirpath}:/input" \
  -v "${refdir}:/ref" \
  -v "${dirpath}:/output" \
  teambraker/braker3:latest \
  braker.pl --genome="/input/${file_name}.masked" --prot_seq="/ref/${ref_name}" --workingdir=/output/braker3 --threads=$THREADS
   --skipOptimize to have it run faster for testing

# Check if the braker.aa file was produced and is not empty
if [ -f "$OUTPUT_DIR/$BASE_NAME/braker3/braker.aa" ] && [ -s "$OUTPUT_DIR/$BASE_NAME/braker3/braker.aa" ]; then
    echo "Braker ran successfully."
else
    echo "Error: braker.aa file not found or is empty in the directory $OUTPUT_DIR/$BASE_NAME/braker3, please recheck your run."
    exit 1
fi

# Copy the braker.aa file to the main directory and rename it
cp "$OUTPUT_DIR/$BASE_NAME/braker3/braker.aa" "$OUTPUT_DIR/$BASE_NAME/${BASE_NAME}_proteome.fa"

# Deactivating the previous conda environment and activating the correct one
conda deactivate
conda activate busco_env
#################################################################################################

# Running BUSCO analysis
cd "$OUTPUT_DIR/$BASE_NAME"
busco -i "$OUTPUT_DIR/$BASE_NAME/${BASE_NAME}_proteome.fa" -o busco_annotated -m protein -l $LINEAGE --cpu $THREADS -f
rm -r "$OUTPUT_DIR/$BASE_NAME/busco_downloads"

# Deactivating the previous conda environment to run interproscan in base
conda deactivate
#################################################################################################

echo "Running Interproscan"
sed -i 's/*//g' "$OUTPUT_DIR/$BASE_NAME/${BASE_NAME}_proteome.fa"
bash ~/my_interproscan/interproscan-5.69-101.0/interproscan.sh -i "$OUTPUT_DIR/$BASE_NAME/${BASE_NAME}_proteome.fa" -f tsv -d "$OUTPUT_DIR/$BASE_NAME" -cpu $THREADS

# Check if the braker.aa file was produced and is not empty
if [ -f "$OUTPUT_DIR/$BASE_NAME/${BASE_NAME}_proteome.fa.tsv" ] && [ -s "$OUTPUT_DIR/$BASE_NAME/${BASE_NAME}_proteome.fa.tsv" ]; then
    echo "Interproscan ran successfully."
else
    echo "Error: proteome.fa.tsv file not found or is empty in the directory $OUTPUT_DIR/$BASE_NAME, please recheck your run."
    exit 1
fi

#################################################################################################

echo "Processing final proteome file"
# Files
input_file="$OUTPUT_DIR/$BASE_NAME/${BASE_NAME}_proteome.fa.tsv"
proteome_file="$OUTPUT_DIR/$BASE_NAME/${BASE_NAME}_proteome.fa"
header_file="$OUTPUT_DIR/$BASE_NAME/${BASE_NAME}_annotated_headers.txt"
output_file="$OUTPUT_DIR/$BASE_NAME/${BASE_NAME}_annotated_proteome.fa"

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
    echo "*               Annotation Pipeline Completed!             *"
    echo "*                                                          *"
    echo "*                        Success!                          *"
    echo "*                                                          *"
    echo "************************************************************"
else
    echo "Error: final annotated fasta file not created or is empty, please recheck your run"
    exit 1
fi

#################################################################################################

# Print the end date and time (append to file)
end_time=$(date "+%Y-%m-%d %H:%M:%S")
end_seconds=$(date +%s)

echo "End Time: $end_time" >> "$log_file"

# Calculating runtime and Converting runtime to human readable
runtime=$((end_seconds - start_seconds))
hours=$((runtime / 3600))
minutes=$(((runtime % 3600) / 60))
seconds=$((runtime % 60))
formatted_runtime=$(printf "%02d:%02d:%02d" $hours $minutes $seconds)

echo "Runtime: $formatted_runtime (HH:MM:SS)" >> "$log_file"

#################################################################################################
