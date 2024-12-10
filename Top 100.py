import pandas as pd
import numpy as np
import os

# Paths for chunk files and SNP files
chunk_paths = [f"/home/hp/Desktop/SNP Analysis/Processed Data/Selected SNP/Output folder 1chunk{i}.csv" for i in range(1, 111)]
snp_paths = [f"/home/hp/Desktop/SNP Analysis/Processed Data/Final_SNP/chunk{i}.csv" for i in range(1, 111)]

# Exclude indices as specified
excluded_indices = {8, 21, 27, 37, 38, 39, 99, 100, 109, 110}
included_indices = [i for i in range(1, 111) if i not in excluded_indices]

# Initialize a matrix (43 rows by length of included indices)
num_rows = 43
mm_matrix = np.empty((num_rows, len(included_indices)), dtype=object)

# Process each valid file
for idx, chunk_idx in enumerate(included_indices):
    # Read the SNP and chunk files
    snp_data = pd.read_csv(snp_paths[chunk_idx - 1], header=None).T  # Transpose SNP data
    chunk_data = pd.read_csv(chunk_paths[chunk_idx - 1])
    
    # Extract position from the first cell in the chunk file
    position_str = chunk_data.iloc[0, 0]
    position_value = int(position_str.split("_")[1])  # Extract numeric position
    
    # Match the position in the SNP file
    position_column = np.where(snp_data.iloc[0, :] == position_value)[0]
    
    # Populate the MM matrix with matching SNP data
    if position_column.size > 0:
        mm_matrix[:, idx] = snp_data.iloc[:, position_column[0]].values
    else:
        print(f"Position {position_value} not found in SNP file for chunk {chunk_idx}")

# Set row names from the SNP file
row_names = snp_data.index.tolist()

# Convert MM matrix into a DataFrame for easier handling
mm_df = pd.DataFrame(mm_matrix, index=row_names)

# Save to a CSV file if needed
mm_df.to_csv("processed_MM_matrix.csv", index=True)
