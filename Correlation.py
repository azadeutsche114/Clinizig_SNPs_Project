import pandas as pd
import numpy as np
import vcf
import matplotlib.pyplot as plt
import seaborn as sns

# Load VCF file
vcf_file = '/media/hp/SVM/1001_GENOME_PROJECT/1001 genomes chunk file/Chunk/vcf_chunk_1.vcf'
vcf_reader = vcf.Reader(filename=vcf_file)

# Extract fixed (meta) data
vcf_meta = pd.DataFrame(vcf_reader.metadata)

# Extract genotype, GQ, DP from VCF file
gt_data = []
gq_data = []
dp_data = []
for record in vcf_reader:
    gt_data.append(record.samples[0]['GT'] if 'GT' in record.samples[0] else None)
    gq_data.append(record.samples[0]['GQ'] if 'GQ' in record.samples[0] else None)
    dp_data.append(record.samples[0]['DP'] if 'DP' in record.samples[0] else None)

# Convert to DataFrame
genotype_data = pd.DataFrame({
    'GT': gt_data,
    'GQ': gq_data,
    'DP': dp_data
})

# Replace NA values
genotype_data['GT'].replace({None: '0'}, inplace=True)
genotype_data['GQ'].fillna(0, inplace=True)
genotype_data['DP'].fillna(0, inplace=True)

# Create combined dataset with fixed meta and genotype data
vcf_df = pd.concat([vcf_meta, genotype_data], axis=1)

# Plotting
# Genotype Frequencies
plt.figure(figsize=(8, 6))
vcf_df['GT'].value_counts().plot(kind='bar', color='skyblue')
plt.title('Genotype Frequencies')
plt.xlabel('Genotype')
plt.ylabel('Frequency')
plt.show()

# Pie chart for Genotype Distribution
plt.figure(figsize=(8, 6))
vcf_df['GT'].value_counts().plot(kind='pie', autopct='%1.1f%%', startangle=90, colors=sns.color_palette('pastel'))
plt.title('Genotype Distribution')
plt.ylabel('')
plt.show()

# Plot Genotype Quality (GQ)
plt.figure(figsize=(8, 6))
sns.histplot(vcf_df['GQ'], kde=True, color='salmon', bins=20)
plt.title('Histogram of Genotype Quality (GQ)')
plt.xlabel('Genotype Quality (GQ)')
plt.ylabel('Count')
plt.show()

# Box plot of Genotype Quality (GQ)
plt.figure(figsize=(8, 6))
sns.boxplot(x=vcf_df['GQ'], color='lightblue')
plt.title('Box Plot of Genotype Quality (GQ)')
plt.ylabel('Genotype Quality (GQ)')
plt.show()

# Plot Read Depth (DP)
plt.figure(figsize=(8, 6))
sns.histplot(vcf_df['DP'], kde=True, color='lightgreen', bins=20)
plt.title('Histogram of Read Depth (DP)')
plt.xlabel('Read Depth (DP)')
plt.ylabel('Count')
plt.show()

# Box plot of Read Depth (DP)
plt.figure(figsize=(8, 6))
sns.boxplot(x=vcf_df['DP'], color='lightcoral')
plt.title('Box Plot of Read Depth (DP)')
plt.ylabel('Read Depth (DP)')
plt.show()

# Correlation heatmaps
# Correlation for GQ data
cor_matrix_gq = vcf_df[['GQ']].corr()
plt.figure(figsize=(8, 6))
sns.heatmap(cor_matrix_gq, annot=True, cmap='coolwarm', fmt='.2f', cbar=True)
plt.title('Correlation Heatmap for Genotype Quality (GQ)')
plt.show()

# Correlation for DP data
cor_matrix_dp = vcf_df[['DP']].corr()
plt.figure(figsize=(8, 6))
sns.heatmap(cor_matrix_dp, annot=True, cmap='coolwarm', fmt='.2f', cbar=True)
plt.title('Correlation Heatmap for Read Depth (DP)')
plt.show()

# Save the DataFrame to CSV
vcf_df.to_csv('/home/hp/Desktop/Chunk 1/chunk1.csv', index=False)
