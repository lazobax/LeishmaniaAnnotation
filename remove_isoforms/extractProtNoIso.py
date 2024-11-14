from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_proteome_no_isoforms(gbff_file, output_fasta):
    print("Processing " + gbff_file)
    proteins = {}
    for record in SeqIO.parse(gbff_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                # Extract the gene name or protein ID
                gene_name = feature.qualifiers.get("locus_tag", ["unknown_gene"])[0]
                protein_id = feature.qualifiers.get("protein_id", [gene_name])[0]
                product = feature.qualifiers.get("product", ["unknown_product"])[0]
                # print(gene_name)
                if "translation" in feature.qualifiers:
                    protein_seq = feature.qualifiers["translation"][0]
                    
                    # If we have already seen this gene or protein ID, compare lengths
                    if gene_name in proteins:
                        existing_protein = proteins[gene_name]
                        if len(protein_seq) > len(existing_protein.seq):
                            print(f"Replacing {gene_name} with longer isoform")
                            proteins[gene_name] = protein_record
                        continue

                    # Create a SeqRecord object
                    protein_record = SeqRecord(
                        Seq(protein_seq),
                        id=protein_id,
                        name=gene_name,
                        description=f"{gene_name} {product}"
                    )
                    proteins[gene_name] = protein_record

    # Write the unique protein sequences to a FASTA file
    SeqIO.write(proteins.values(), output_fasta, "fasta")


# Read species names from filenames.txt
with open("filenames.txt", "r") as file:
    species_list = file.read().splitlines()

for sp in species_list:
    # Specify the input GBFF file and output FASTA file
    gbff_file = f"gbffs/{sp}.gbff"
    output_fasta = f"{sp}noiso.fasta"

    # Extract the proteome without isoforms
    extract_proteome_no_isoforms(gbff_file, output_fasta)
