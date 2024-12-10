import pandas as pd
import numpy as np
import os
import vcfpy

# Load accession data from Excel file
accessions = pd.read_excel("Desktop/SNP Analysis/Accessions shortlisted/Assessions.xlsx", header=None, skiprows=1)
access_no = accessions.iloc[0, :].dropna().unique()

# Load supplementary table
supplementary_table = pd.read_csv("~/Desktop/SNP Analysis/Accessions shortlisted/Supplementary Table S1.csv")
w1 = supplementary_table[supplementary_table["Accession_Code"].isin(access_no)].index

# Paths for input and output folders
vcf_folder = "/media/hp/SVM/1001_GENOME_PROJECT/1001 genomes chunk file/Chunk"
output_folder1 = "/home/hp/Desktop/SNP Analysis/Processed Data/GQ/"
output_folder2 = "/home/hp/Desktop/SNP Analysis/Processed Data/DP/"

# Get a list of all VCF files in the folder
vcf_files = [os.path.join(vcf_folder, f) for f in os.listdir(vcf_folder) if f.endswith(".vcf")]

# Process each VCF file
for i, vcf_file in enumerate(vcf_files, start=1):
    print(f"Processing file {i}: {vcf_file}")
    
    # Read VCF file
    reader = vcfpy.Reader.from_path(vcf_file)
    
    # Initialize GQ and DP data dictionaries
    gq_data = {}
    dp_data = {}
    
    # Process records in the VCF file
    for record in reader:
        for sample in record.call_for_sample.keys():
            call = record.call_for_sample[sample]
            # Extract GQ and DP
            gq_data.setdefault(sample, []).append(call.data.get("GQ", 0))
            dp_data.setdefault(sample, []).append(call.data.get("DP", 0))
    
    # Convert GQ and DP data to DataFrames
    gq_df = pd.DataFrame(gq_data)
    dp_df = pd.DataFrame(dp_data)
    
    # Match accession codes
    w3 = supplementary_table.loc[w1, "Accession_Code"].astype(str).map(lambda x: x in gq_df.columns)
    gq_df = gq_df.loc[:, w3]
    dp_df = dp_df.loc[:, w3]
    
    # Replace NaN with 0
    gq_df.fillna(0, inplace=True)
    dp_df.fillna(0, inplace=True)
    
    # Save GQ and DP data to CSV files
    gq_output_filename1 = os.path.join(output_folder1, f"chunk{i}.csv")
    dp_output_filename2 = os.path.join(output_folder2, f"chunk{i}.csv")
    gq_df.to_csv(gq_output_filename1, index=False)
    dp_df.to_csv(dp_output_filename2, index=False)
    
    print(f"Processed chunk {i}")
