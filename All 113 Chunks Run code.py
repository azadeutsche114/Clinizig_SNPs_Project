import os
import pandas as pd
import numpy as np
import vcfpy

# Load Excel and CSV files
accessions = pd.read_excel("Desktop/SNP Analysis/Accessions shortlisted/Assessions.xlsx", header=None, skiprows=1)
access_no = accessions.iloc[0].dropna().unique().astype(int)

st = pd.read_csv("~/Desktop/SNP Analysis/Accessions shortlisted/Supplementary Table S1.csv")
w1 = st[st['Accession_Code'].isin(access_no)].index

# Path to the folder containing VCF files
vcf_folder = "/media/hp/SVM/1001_GENOME_PROJECT/1001 genomes chunk file/Chunk"
vcf_files = [os.path.join(vcf_folder, f) for f in os.listdir(vcf_folder) if f.endswith(".vcf")]

# Initialize an empty list to store results
results_list = []

# Function to extract GQ and DP data from VCF
def extract_gt_data(vcf_reader, field, samples):
    data = []
    for record in vcf_reader:
        sample_values = [
            record.call_for_sample[sample].data.get(field) for sample in samples
        ]
        data.append(sample_values)
    return np.array(data, dtype=float)

# Loop through each VCF file
for vcf_file in vcf_files:
    with vcfpy.Reader.from_path(vcf_file) as reader:
        samples = reader.header.samples.names
        # Match accessions
        matched_samples = [samples[i] for i in range(len(samples)) if int(samples[i]) in access_no]
        
        # Extract GQ and DP data
        gq_data = extract_gt_data(reader, "GQ", matched_samples)
        dp_data = extract_gt_data(reader, "DP", matched_samples)
        
        # Replace NaN with 0
        gq_data = np.nan_to_num(gq_data, nan=0.0)
        dp_data = np.nan_to_num(dp_data, nan=0.0)
        
        # Calculate mean and standard deviation
        gq_mean = np.mean(gq_data, axis=0)
        gq_std = np.std(gq_data, axis=0)
        dp_mean = np.mean(dp_data, axis=0)
        dp_std = np.std(dp_data, axis=0)
        
        # Store results
        results_df = pd.DataFrame({
            'File': os.path.basename(vcf_file),
            'Accession_Code': matched_samples,
            'GQ_mn': gq_mean,
            'GQ_sd': gq_std,
            'DP_mn': dp_mean,
            'DP_sd': dp_std
        })
        results_list.append(results_df)

# Combine all results into a single DataFrame
final_results = pd.concat(results_list, ignore_index=True)

# Save the results to a CSV file
final_results.to_csv("GQ_DP_statistics.csv", index=False)

# Further processing with GQ_DP_statistics
gq_dp_statistics = pd.read_csv("GQ_DP_statistics.csv")
unique_accessions = gq_dp_statistics['Accession_Code'].unique()
gq_dp_final = pd.DataFrame(index=unique_accessions, columns=['GQ_mn', 'GQ_sd', 'DP_mn', 'DP_sd'])

for accession in unique_accessions:
    accession_data = gq_dp_statistics[gq_dp_statistics['Accession_Code'] == accession]
    gq_dp_final.loc[accession] = accession_data.iloc[:, 2:6].mean().values

# Adding additional columns from the Supplementary Table
gq_dp_final['R.1'] = st.loc[w1, 'R.1'].values
gq_dp_final['R.2'] = st.loc[w1, 'R.2'].values
gq_dp_final['R.3'] = st.loc[w1, 'R.3'].values

# Save final results to a CSV file
gq_dp_final.to_csv("GQ_DP_final.csv", index=False)
