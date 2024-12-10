import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
from flask import Flask, render_template_string, send_file

# Read Excel and CSV data
assessions = pd.read_excel("Desktop/SNP Analysis/Accessions shortlisted/Assessions.xlsx", header=None, skiprows=1)
access_no = pd.unique(assessions.iloc[0, :])

st = pd.read_csv("~/Desktop/SNP Analysis/Accessions shortlisted/Supplementary Table S1.csv")
w1 = st[st['Accession_Code'].isin(access_no)].index

# Define the path to the directory containing VCF files
vcf_dir = "/media/hp/SVM/1001_GENOME_PROJECT/1001 genomes chunk file/Chunk/"

# List all VCF files in the directory
vcf_files = glob.glob(os.path.join(vcf_dir, "*.vcf"))

# Initialize a dictionary to store unique values from each file
all_unique_values = {}
iteration = 1

# Iterate over each VCF file
for vcf_file in vcf_files:
    print(f"Processing Iteration: {iteration}")
    # Read VCF file (assuming it has genotype data in tabular format)
    vcf_data = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
    
    # Extract genotype data (assuming GT column exists)
    gt_data = vcf_data.iloc[:, -1]  # Replace with actual genotype column index if known
    
    # Map selected accession codes to columns (modify as needed for real data)
    gt_data_selected = gt_data.loc[:, st.loc[w1, 'Accession_Code']]
    
    # Convert to DataFrame
    gt_data_df = pd.DataFrame(gt_data_selected)
    
    # Get unique values
    unique_values = gt_data_df.stack().value_counts().to_dict()
    
    # Store unique values in the dictionary
    all_unique_values[vcf_file] = unique_values
    iteration += 1

# Generate and save 112 full-sized bar plots
os.makedirs("plots", exist_ok=True)
for i, (file, values) in enumerate(all_unique_values.items(), start=1):
    plot_data = pd.DataFrame(list(values.items()), columns=["Genotype", "Frequency"])
    
    plt.figure(figsize=(5, 4))
    plt.bar(plot_data["Genotype"], plot_data["Frequency"], color='blue')
    plt.title(f"Plot {i}")
    plt.xlabel("Genotype")
    plt.ylabel("Frequency")
    plt.savefig(f"plots/plot_{i}.png", dpi=100)
    plt.close()

# Flask App for Interactive Plot Display
app = Flask(__name__)

# HTML template for the Flask app
template = """
<!DOCTYPE html>
<html>
<head>
    <title>Interactive Plot Selector</title>
</head>
<body>
    <h1>Interactive Plot Selector</h1>
    <div>
        {% for i in range(1, 113) %}
            <a href="/plot/{{ i }}">Plot {{ i }}</a><br>
        {% endfor %}
    </div>
    <div>
        {% if selected_plot %}
            <h2>Selected Plot: {{ selected_plot }}</h2>
            <img src="/static/plots/plot_{{ selected_plot }}.png" alt="Plot {{ selected_plot }}">
        {% endif %}
    </div>
</body>
</html>
"""

@app.route("/")
@app.route("/plot/<int:plot_id>")
def index(plot_id=None):
    return render_template_string(template, selected_plot=plot_id)

if __name__ == "__main__":
    app.run(debug=True)
