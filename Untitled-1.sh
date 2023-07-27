#!/bin/bash

# Path to the input fasta file containing peptides
peptides_file="path/to/your/peptides.fasta"

# Create a temporary directory to store filtered fasta files
tmp_dir=$(mktemp -d)

# Use seqkit to split the fasta file into separate files based on peptide length
seqkit split -O "$tmp_dir" -i -s "^([ACTGactg]+)$" "$peptides_file"

# Loop through each filtered fasta file and run pvacbind
for filtered_file in "$tmp_dir"/*.fasta; do
    peptide_length=$(basename "$filtered_file" .fasta)
    pvacbind run --iedb-install-directory /opt/iedb "$filtered_file" "$meta.ID" "$Class_I" -e1 "$peptide_length" "NetMHCpan"
done

# Clean up temporary directory
rm -rf "$tmp_dir"

pvacbind run --iedb-install-directory /opt/iedb $peptides_I $meta.ID \$Class_I -e1 9 "NetMHCpan" . 
