#--------------------------------------------------------------
#
# ILC Jun 2025
# TFM - Universidad Internacional de Valencia - Máster en Bioinformática
#
#--------------------------------------------------------------
#
# GWAS PLOT CREATION
#
# -------------------------------------------------------------
#

import pandas as pd
import os
import re
import glob
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#--------------------------------------------------------------
#
# Required input and output files.

output_dir = "07_plots/"
input_files = glob.glob("04_annotation_results/result_ph_*_annotated.tsv")
adapted_files = glob.glob("07_plots/*_gwaslab.tsv")
log_file = "output_std_07.txt"

# Create directory

os.makedirs(output_dir, exist_ok=True)

# Redirect std and err output to log file.

sys.stdout = open (log_file, 'w')
sys.stderr = sys.stdout

#
#--------------------------------------------------------------
#
# 1. File adaptation for GWASlab.
#

def add_unique_snp_column(input_file, output_dir):
    print(f"Processing {input_file} file.")
    df = pd.read_csv(input_file, sep='\t')

    df['UNIQSNP'] = df['CHR'].astype(str) + ":" + df['BP'].astype(str) + "_" + df['A1']

    columns_order = ['UNIQSNP', 'CHR', 'BP', 'A1', 'P', 'BETA']
    df = df[[col for col in columns_order if col in df.columns]]

    match = re.search(r"result_ph_(.+?)_annotated\.tsv", os.path.basename(input_file))
    if not match:
        print(f"Error: could not extract phenotype from {input_file} file.")
        return None

    phenotype = match.group(1)
    output_file = os.path.join(output_dir, f"{phenotype}_gwaslab.tsv")
    
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved adapted file as {output_file}")

    return output_file

#
#--------------------------------------------------------------
#
# 2. Manhhatan plot.
#

def generate_manhattan (adapted_file):
    phenotype = re.search(r"(.+?)_gwaslab\.tsv", os.path.basename(adapted_file)).group(1)
    output_file = os.path.join(output_dir, f"{phenotype}_manhattanplot.png")
    title = f"Manhhatan Plot - {phenotype}"

    print(f"Reading {adapted_file}")
    print(f"Saved in {output_file}")

    try:
        df = pd.read_csv(adapted_file, sep='\t')
        df['CHR'] = df['CHR'].astype(int)
        df['BP'] = df['BP'].astype(int)
        df['P'] = pd.to_numeric(df['P'], errors= 'coerce')
        df = df.dropna(subset=['P'])

        df['-log10(P)'] = -np.log10(df['P'])

        df = df.sort_values(['CHR', 'BP'])
        df['ind'] = range(len(df))

        colors = ['#7f7f7f', '#ff7f0e']
        df['color'] = df['CHR'] % 2
        df['color'] = df['color'].apply(lambda x: colors[x])
       
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.scatter(df['ind'], df['-log10(P)'], c=df['color'], s=10, alpha=0.6, edgecolor='none')

        ax.axhline(-np.log10(5e-8), color='grey', linestyle='--', lw=1)

        ticks = df.groupby('CHR')['ind'].median()
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks.index)

        ax.set_xlabel('Chromosome')
        ax.set_ylabel('-Log10(P)')
        ax.set_title(title)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()

        print(f"Finished: {phenotype}\n")

    except Exception as e: 
        print(f"Error in processing {adapted_file}: {e}\n")

#
#--------------------------------------------------------------
#
# 3. QQ plot.
#

def generate_qq (adapted_file):
    phenotype = re.search(r"(.+?)_gwaslab\.tsv", os.path.basename(adapted_file)).group(1)
    output_file = os.path.join(output_dir, f"{phenotype}_QQplot.png")
    title = f"QQ Plot - {phenotype}"

    print(f"Reading {adapted_file}")
    print(f"Saved in {output_file}")

    try:
        df = pd.read_csv(adapted_file, sep='\t')

        df['P'] = pd.to_numeric(df['P'], errors='coerce')
        df = df.dropna(subset=['P'])
        df = df[df['P'] > 0]  

        n = len(df)
        expected = -np.log10(np.linspace(1 / n, 1, n))
        observed = -np.log10(np.sort(df['P']))

        plt.figure(figsize=(6, 6))
        plt.plot(expected, observed, 'o', markersize=2, alpha=0.6,color='#ff7f0e', label='Observed vs. Expected')
        plt.plot([0, max(expected)], [0, max(expected)], color='#7f7f7f', linestyle='--', label='Expected')
        plt.xlabel('Expected -log10(P)')
        plt.ylabel('Observed -log10(P)')
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()

        print(f"Finished: {phenotype}\n")

    except Exception as e:
        print(f"Error in processing {adapted_file}: {e}\n")

#
#-----------------------------------------------------------
#

for input_file in input_files:
    adapted = add_unique_snp_column(input_file, output_dir)
    if adapted: 
        generate_manhattan(adapted)
        generate_qq(adapted)

print("All plots successfully generated.")

#
#-----------------------------------------------------------
#

