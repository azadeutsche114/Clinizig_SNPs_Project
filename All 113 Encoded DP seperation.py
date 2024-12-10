import pandas as pd
import os

# Define input and output folders
gq_folder = "/home/hp/Desktop/SNP Analysis/Processed Data/GQ"
dp_folder = "/home/hp/Desktop/SNP Analysis/Processed Data/Encoded_DP/"
output_folder = "/home/hp/Desktop/SNP Analysis/Processed Data/Final DP_GQ Encoded/"

# Get a list of files
gq_files = sorted([os.path.join(gq_folder, f) for f in os.listdir(gq_folder) if f.endswith('.csv')])
dp_files = sorted([os.path.join(dp_folder, f) for f in os.listdir(dp_folder) if f.endswith('.csv')])

# Process each pair of GQ and DP files
for i, (gq_file, dp_file) in enumerate(zip(gq_files, dp_files), start=1):
    # Read the CSV files
    gq_df = pd.read_csv(gq_file, header=None)
    dp_df = pd.read_csv(dp_file, header=None)
    
    # Assign column headers from the first row and remove it
    gq_df.columns = gq_df.iloc[0]
    dp_df.columns = dp_df.iloc[0]
    gq_df = gq_df[1:]
    dp_df = dp_df[1:]
    
    # Convert values to numeric
    gq_df = gq_df.apply(pd.to_numeric)
    dp_df = dp_df.apply(pd.to_numeric)
    
    # Perform element-wise multiplication
    multiplied_df = gq_df * dp_df
    
    # Define the output file path
    output_file = os.path.join(output_folder, f"chunk{i}.csv")
    
    # Save the resulting data to a new file
    multiplied_df.to_csv(output_file, index=False, header=False)
    
    print(f"Processed chunk {i}")  # Status update
