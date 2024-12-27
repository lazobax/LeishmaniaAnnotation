#!/bin/bash
set -euo pipefail

################################################################################
# run_annotation_pipeline.sh
#
# Orchestrates:
#   1. Downloading genome/protein data from NCBI (via datasets)
#   2. Removing isoforms with a Python script
#   3. Running InterProScan for annotation
#
# Usage: ./run_annotation_pipeline.sh <csv_file> <threads> <output_dir>
#
#   <csv_file>   CSV with columns: accession, species_name (no header)
#   <threads>    Number of CPU threads for InterProScan
#   <output_dir> Directory to store downloads and final annotation results
################################################################################

if [[ "$#" -lt 3 ]]; then
    echo "Usage: $0 <csv_file> <threads> <output_dir>"
    exit 1
fi

CSV_FILE="$1"
THREADS="$2"
OUTPUT_DIR="$3"

# Create the main output directory if it doesn’t exist
mkdir -p "${OUTPUT_DIR}"

# Log file for any download/rename errors
ERROR_LOG="${OUTPUT_DIR}/download_errors.txt"
rm -f "${ERROR_LOG}"

################################################################################
# 1. Define helper functions for NCBI download + rename
################################################################################

download_genome() {
    local species="$1"
    local accession="$2"

    # Clean species name for filesystem
    local safe_species
    safe_species=$(echo "$species" | tr -d '\n' | xargs | sed 's/[^a-zA-Z0-9]/_/g')

    echo "Downloading for species: $safe_species and accession: $accession"

    # Prepare a local directory to hold the data temporarily
    local species_tmp_dir="${OUTPUT_DIR}/${safe_species}_dataset"
    mkdir -p "${species_tmp_dir}"

    local zip_file="${OUTPUT_DIR}/${safe_species}.zip"

    if ! datasets download genome accession "$accession" --include "gbff,protein" --filename "${zip_file}"; then
        echo "[ERROR] $safe_species: Failed to download dataset for accession $accession" >> "${ERROR_LOG}"
        return 1
    fi

    # Unzip
    if ! unzip -q "${zip_file}" -d "${species_tmp_dir}"; then
        echo "[ERROR] $safe_species: Failed to unzip dataset" >> "${ERROR_LOG}"
        return 1
    fi

    return 0
}

rename_files() {
    local species="$1"
    local accession="$2"

    local safe_species
    safe_species=$(echo "$species" | tr -d '\n' | xargs | sed 's/[^a-zA-Z0-9]/_/g')

    local species_tmp_dir="${OUTPUT_DIR}/${safe_species}_dataset"
    local accession_path="${species_tmp_dir}/ncbi_dataset/data/${accession}"

    local gbffFile="${accession_path}/genomic.gbff"
    local protFile="${accession_path}/protein.faa"

    # Destination subfolders
    local gbff_dir="${OUTPUT_DIR}/gbff"
    local proteomes_dir="${OUTPUT_DIR}/proteomes"
    mkdir -p "${gbff_dir}" "${proteomes_dir}"

    # Move files to final location
    if ! mv "${gbffFile}" "${gbff_dir}/${safe_species}.gbff"; then
        echo "[ERROR] $safe_species: Failed to move gbff file" >> "${ERROR_LOG}"
        return 1
    fi

    if ! mv "${protFile}" "${proteomes_dir}/${safe_species}.fasta"; then
        echo "[ERROR] $safe_species: Failed to move protein file" >> "${ERROR_LOG}"
        return 1
    fi

    # Cleanup
    rm -rf "${species_tmp_dir}"
    rm -f "${OUTPUT_DIR}/${safe_species}.zip"

    return 0
}

################################################################################
# 2. Create the Python script to remove isoforms
################################################################################

cat <<EOF > removeIsoforms.py
#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_proteome_no_isoforms(gbff_file, output_fasta):
    """
    Keeps only the longest isoform per gene from a GenBank (.gbff) file
    and writes these sequences to the specified FASTA file.
    """
    print(f"Processing {gbff_file}")
    proteins = {}

    for record in SeqIO.parse(gbff_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                # Use gene name (if missing, fallback to "unknown_gene")
                gene_name = feature.qualifiers.get("locus_tag", ["unknown_gene"])[0]
                # Use protein_id if present, else fallback to gene_name
                protein_id = feature.qualifiers.get("protein_id", [gene_name])[0]
                # Use the product description if present, else fallback to "unknown_product"
                product = feature.qualifiers.get("product", ["unknown_product"])[0]

                # Check if there's an actual protein translation
                if "translation" in feature.qualifiers:
                    protein_seq = feature.qualifiers["translation"][0]
                    
                    # If we've seen this gene before, see if new isoform is longer
                    if gene_name in proteins:
                        existing_protein = proteins[gene_name]
                        if len(protein_seq) > len(existing_protein.seq):
                            print(f"Replacing isoform for gene '{gene_name}' "
                                  f"({len(existing_protein.seq)} aa) with longer isoform "
                                  f"({len(protein_seq)} aa).")
                            protein_record = SeqRecord(
                                Seq(protein_seq),
                                id=protein_id,
                                name=gene_name,
                                description=f"{gene_name} {product}"
                            )
                            proteins[gene_name] = protein_record
                    else:
                        # First isoform for this gene—store it
                        protein_record = SeqRecord(
                            Seq(protein_seq),
                            id=protein_id,
                            name=gene_name,
                            description=f"{gene_name} {product}"
                        )
                        proteins[gene_name] = protein_record

    # Write one protein per gene (the longest isoform) to the FASTA output
    SeqIO.write(proteins.values(), output_fasta, "fasta")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: removeIsoforms.py <input_gbff> <output_fasta>")
        sys.exit(1)

    gbff_file = sys.argv[1]
    output_fasta = sys.argv[2]

    extract_proteome_no_isoforms(gbff_file, output_fasta)
EOF

chmod +x removeIsoforms.py

################################################################################
# 3. Define a function to run InterProScan, adapting the provided script
################################################################################

run_interproscan() {
    local species="$1"
    local threads="$2"
    local output_dir="$3"

    local safe_species
    safe_species=$(echo "$species" | tr -d '\n' | xargs | sed 's/[^a-zA-Z0-9]/_/g')

    # Where is the no-isoform FASTA located?
    local noiso_fasta="${output_dir}/proteomes/${safe_species}_noiso.fasta"
    local species_outdir="${output_dir}/annotated_proteomes/${safe_species}"
    mkdir -p "${species_outdir}"

    echo "Running InterProScan for ${safe_species}"

    # REPLACE !!! BELOW WITH /home/<user_name> BECAUSE FOR SOME REASON ~ DIDN'T WORK
    # Example path to your InterProScan installation:
    local interproscan_bin="/home/lazo_ozal/my_interproscan/interproscan-5.69-101.0/interproscan.sh"

    bash "${interproscan_bin}" \
        -i "${noiso_fasta}" \
        -f tsv \
        -d "${species_outdir}" \
        -cpu "${threads}"

    local tsv_file="${species_outdir}/${safe_species}_noiso.fasta.tsv"
    if [[ -f "${tsv_file}" && -s "${tsv_file}" ]]; then
        echo "InterProScan completed successfully for ${safe_species}."
    else
        echo "[ERROR] No valid TSV output found for ${safe_species} in ${species_outdir}"
        exit 1
    fi

    ############################################################################
    # Process final proteome file
    ############################################################################

    local input_file="${tsv_file}"
    local proteome_file="${noiso_fasta}"
    local header_file="${species_outdir}/${safe_species}_annotated_headers.txt"
    local output_file="${species_outdir}/${safe_species}_annotated_proteome.fa"

    echo "Post-processing InterProScan annotations for ${safe_species}"

    # 1) Extract best annotation lines, storing them in header_file
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

        # Keep the line with the lowest e-value for each ID
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
    }' "${input_file}" | sort -k1,1 > "${header_file}"

    # 2) For any IDs in the proteome not in the TSV, add them to the header file
    comm -13 \
        <(awk -F'\t' '{print ">" $1}' "${input_file}" | sort -u) \
        <(grep ">" "${proteome_file}" | sort -u) >> "${header_file}"

    # 3) Create a temporary file to map ID -> annotated header
    local temp_file
    temp_file=$(mktemp)
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
    }' "${header_file}" > "${temp_file}"

    # 4) Replace headers in the no-isoform proteome with the best annotation
    awk -v header_file="${temp_file}" '
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
    }' "${proteome_file}" > "${output_file}"

    rm -f "${temp_file}"

    # Final check
    if [[ -s "${output_file}" ]]; then
        echo "************************************************************"
        echo "* Annotation for ${safe_species} completed successfully!   *"
        echo "* Final proteome: ${output_file}"
        echo "************************************************************"
    else
        echo "[ERROR] No final annotated proteome for ${safe_species}"
        exit 1
    fi
}

################################################################################
# 4. Main loop to process each line in the CSV
################################################################################

while IFS=',' read -r accession species; do
    # Skip empty or comment lines
    [[ -z "${accession}" || -z "${species}" ]] && continue

    echo "=================================================================="
    echo "PROCESSING: ${species} (Accession: ${accession})"
    echo "=================================================================="

    # 1. Download
    if ! download_genome "${species}" "${accession}"; then
        echo "Skipping ${species} due to download error"
        continue
    fi

    # 2. Rename
    if ! rename_files "${species}" "${accession}"; then
        echo "Skipping ${species} due to rename error"
        continue
    fi

    # 3. Remove isoforms
    #    - The newly placed gbff file is at:  <output_dir>/gbff/<species_safename>.gbff
    #    - We will produce a no-isoform FASTA at: <output_dir>/proteomes/<species_safename>_noiso.fasta

    safe_species=$(echo "$species" | tr -d '\n' | xargs | sed 's/[^a-zA-Z0-9]/_/g')
    gbff_file="${OUTPUT_DIR}/gbff/${safe_species}.gbff"
    noiso_fasta="${OUTPUT_DIR}/proteomes/${safe_species}_noiso.fasta"

    echo "Removing isoforms for ${species}"
    ./removeIsoforms.py "${gbff_file}" "${noiso_fasta}"

    # 4. Run InterProScan (annotation)
    run_interproscan "${species}" "${THREADS}" "${OUTPUT_DIR}"

    echo "Finished processing ${species}"
    echo
done < "${CSV_FILE}"

echo "All species processed. Check '${OUTPUT_DIR}/annotated_proteomes/' for final results."
