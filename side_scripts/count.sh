#!/bin/bash
# Counts the percentage of hypothetical proteins in a fasta file

# Check if the script is being run in a directory with .fasta files
if ! ls *.fasta 1> /dev/null 2>&1; then
    echo "No .fasta files found in the current directory."
    exit 1
fi

# Loop through each .fasta file in the current directory
for file in *.fasta; do
    # Count occurrences of "hypothetical" (case insensitive)
    hypothetical_count=$(grep -i "hypothetical" "$file" | wc -l)
    
    # Count the number of header lines (lines starting with ">")
    header_count=$(grep -c "^>" "$file")
    
    # Calculate the percentage if headers are found
    if [ "$header_count" -gt 0 ]; then
        percentage=$(echo "scale=2; $hypothetical_count / $header_count * 100" | bc)
        echo "$file: $hypothetical_count hypothetical(s) out of $header_count header(s) = $percentage%"
    else
        echo "$file: No headers found"
    fi
done
