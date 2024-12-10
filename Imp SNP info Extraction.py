import pandas as pd
import os
import vcfpy

# Define the path to your folder containing the VCF files
vcf_folder = "/media/hp/SVM/1001_GENOME_PROJECT/1001 genomes chunk file/Chunk"
output_folder = "/home/hp/Desktop/SNP Analysis/Processed Data/Final_SNP/"

# Create a list of file paths for chunks
chunk_paths = [f"/home/hp/Desktop/SNP Analysis/Processed Data/Selected SNP/Output folder 1chunk{i}.csv" for i in range(71, 86)]

# Get a list of all VCF files in the folder
vcf_files = [os.path.join(vcf_folder, file) for file in os.listdir(vcf_folder) if file.endswith('.vcf')]

# Indices to exclude
excluded_indices = {8, 21, 27, 37, 38, 39, 99, 100, 109, 110}
included_indices = [i for i in range(1, 111) if i not in excluded_indices]

# Loop through each VCF file
for i in range(len(included_indices)):
    # Read the corresponding VCF file
    vcf_file_path = vcf_files[included_indices[i] - 1]  # Adjust index for zero-based Python
    reader = vcfpy.Reader.from_path(vcf_file_path)

    # Load the corresponding chunk file
    chunk_path = chunk_paths[i]
    chunk_data = pd.read_csv(chunk_path)

    # Extract the position value from the first row
    pos_str = chunk_data.iloc[0, 0]
    pos_value = int(pos_str.split("_")[1])  # Extract numeric position

    # Search for the SNP with the given position in the VCF file
    snp_found = None
    for record in reader:
        if record.POS == pos_value:
            snp_found = record
            break

    # Print the metadata of the found SNP
    if snp_found:
        print(f"SNP found in file {vcf_file_path} at position {pos_value}:")
        print(snp_found)
    else:
        print(f"No SNP found in file {vcf_file_path} at position {pos_value}")
