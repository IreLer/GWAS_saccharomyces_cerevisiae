#--------------------------------------------------------------
#
# ILC Jun 2025
# TFM - Universidad Internacional de Valencia - Máster en Bioinformática
#
#--------------------------------------------------------------
#
# VARIANT FILTERING
# 35 PLINK results files.
#
#-------------------------------------------------------------
#
import os
import pandas as pd
import sys
import glob
#
#--------------------------------------------------------------
#
# Required input files.

variant_matrix = "01_modified_files/GWAS_Matrix.bim"
result_files = glob.glob("./02_GWAS_results/result_ph_*.assoc.linear")
output_dir = "03_filtering_results"
log_file = "output_std_03.txt"

# Create directory.

os.makedirs(output_dir, exist_ok=True)

# Redirect std and err output to log file.

sys.stdout = open (log_file, 'w') 
sys.stderr = sys.stdout

#
#--------------------------------------------------------------
#
# 1. Bonferroni correction threshold calculation. 
#

with open(variant_matrix, 'r') as f: 
    n_test = sum(1 for _ in f)

pvalue_threshold = 0.05/n_test

print(f"Number of test: {n_test}")
print(f"Bonferroni-corrected p-value threshold: {pvalue_threshold}.")

#
#--------------------------------------------------------------
#
# 2. Filtering PLINK results by p-value.
#

for file in result_files:
    print(f"Processing {file}...")
    df = pd.read_csv(file, sep=r'\s+')            #Tab separation.
    df_filtered = df[df['P'] < pvalue_threshold]  # Filter by p-value.
    base_name = os.path.basename(file).replace('.assoc.linear','_filtered.assoc.linear')
    output_file = os.path.join(output_dir, base_name)    
    df_filtered.to_csv(output_file, sep='\t', index=False)
    print(f"Filtered file saved as: {output_file}.")

filtered_files = glob.glob("./03_filtering_results/result_ph_*_filtered.assoc.linear")

for file in filtered_files: 
    with open(file, 'r') as f: 
        num_lines = sum(1 for _ in f)

print(f"The file {file} contains {num_lines} lines.")

#
#--------------------------------------------------------------
#
