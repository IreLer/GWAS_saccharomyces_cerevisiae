#--------------------------------------------------------------
#
# ILC Jun 2025
# TFM - Universidad Internacional de Valencia - Máster en Bioinformática
#
#--------------------------------------------------------------
#
# RESULT ANALYSIS
#
# ---------------
#
import pandas as pd
import glob
import re
import os
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#--------------------------------------------------------------
#
# Required input and output files.

output_dir = "05_analysis_results"
annot_files = glob.glob("./04_annotation_results/result_ph_*_annotated.tsv")
betapos_files = glob.glob(os.path.join(output_dir, "*_BETApos.tsv"))
betaneg_files = glob.glob(os.path.join(output_dir, "*_BETAneg.tsv"))
matrix_files = glob.glob(os.path.join(output_dir, "*_matrix_BETA*.tsv"))
log_file = "output_std_05.txt"

# Create directory

os.makedirs(output_dir, exist_ok=True)

# Redirect std and err output to log file.

sys.stdout = open (log_file, 'w')
sys.stderr = sys.stdout

#
#--------------------------------------------------------------
#
# 1. Split BETA positive and negative genes.
#

print(f"Starting BETA preparation...")

for file in annot_files:
    match = re.search(r"result_ph_(.+?)_annotated\.tsv", file)
    
    phenotype = match.group(1)
   
    print(f"Processing phenotype: {phenotype}")

    try:
        df = pd.read_csv(file, sep="\t")
    except Exception as e:
        print(f"Error reading {file}:{e}")
        continue

    df = df[pd.to_numeric(df["BETA"], errors="coerce").notnull()]
    df["BETA"] = df["BETA"].astype(float)

    df_pos = df[df["BETA"]>0]
    df_neg = df[df["BETA"]<0]

    pos_path = os.path.join(output_dir, f"{phenotype}_BETApos.tsv")
    neg_path = os.path.join(output_dir, f"{phenotype}_BETAneg.tsv")
    
    df_pos.to_csv(pos_path, sep="\t", index=False)
    df_neg.to_csv(neg_path, sep="\t", index=False)

    print(f"Written {len(df_pos)} rows to {pos_path}")
    print(f"Written {len(df_neg)} rows to {neg_path}")

print(f"All BETA files have been susccesfully generated.")


#--------------------------------------------------------------------------
#
# 2. Create the genotype - phenotype matrix
#        

ph_1 = ["YPD14", "YPD40", "YPD42"]
ph_2 = ["YPACETATE", "YPETHANOL", "YPGALACTOSE", "YPGLYCEROL", "YPRIBOSE", "YPSORBITOL", "YPXYLOSE"]
ph_3 = ["YPDKCL2M", "YPDNACL15M", "YPDNACL1M"]
ph_4 = ["YPDCUSO410MM", "YPDLICL250MM", "YPDMV","YPDSODIUMMETAARSENITE"]
ph_5 = ["YPDANISO10", "YPDANISO20", "YPDANISO50", "YPDCHX05", "YPDCHX1", "YPD6AU", "YPDHU"]
ph_6 = ["YPDBENOMYL200", "YPDBENOMYL500", "YPDFLUCONAZOLE", "YPDNYSTATIN","YPDCAFEIN40", "YPDCAFEIN50", "YPDDMSO", "YPDETOH", "YPDFORMAMIDE4", "YPDFORMAMIDE5", "YPSDS"]

phenotype_groups = {"heat_stress" : ph_1,
                    "carbon_replacement" : ph_2,
                    "osmotic_stress" : ph_3,
                    "oxidative_stress" : ph_4,
                    "synthesis_inhibitors" : ph_5,
                    "chemical_agents" : ph_6}
                          
def build_matrix(filtered_files, output_filename, phenotypes_of_interest):
    matrix = {}
    all_genes = set()

    for file in filtered_files: 
        match = re.search(r"/([A-Z0-9]+)_BETA(pos|neg)\.tsv$", file)
        if not match: 
            print(f"Skipping file: {file} - not found.")
            continue

        phenotype = match.group(1)
        
        print(f"Detected phenotype: {phenotype}")

        if phenotype not in phenotypes_of_interest: 
            print(f"Ignoring phenotype {phenotype}")
            continue

        try: 
            df =  pd.read_csv(file, sep="\t", usecols=["gene_name", "P"])
        except Exception as e: 
            print(f"Error reading {file}:{e}")
            continue

        df = df[df["P"]>0]
        summary = df.groupby("gene_name")["P"].min()
        matrix[phenotype] = summary
        all_genes.update(summary.index)

    genes_sorted = sorted(all_genes)
    df_matrix = pd.DataFrame(index=genes_sorted)

    for pheno in phenotypes_of_interest:
        data = matrix.get(pheno, {})
        column = [data.get(g,0) for g in df_matrix.index]
        df_matrix[pheno] = column
    
    df_matrix.replace(0,1, inplace=True)
    df_matrix.reset_index(inplace=True)
    df_matrix.rename(columns={"index":"gene"}, inplace=True)

    output_path = os.path.join(output_dir, output_filename)
    df_matrix.to_csv(output_path, sep="\t", index=False)
    print(f"{output_path} matrix generated.")

for group_name, phenotypes in phenotype_groups.items():
    build_matrix(betapos_files, f"{group_name}_matrix_BETApos.tsv", phenotypes)
    build_matrix(betaneg_files, f"{group_name}_matrix_BETAneg.tsv", phenotypes)

print("Matrix susccesfully generated.")

#
#--------------------------------------------------------------------------
#
# 3. Generate the HEATMAP.
#

def generate_heatmap (matrix_path, heatmap_path, title=None):
    try:
        df_matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)

        df_log = -np.log10(df_matrix)

        df_log.replace([np.inf, -np.inf], 300, inplace=True)

        df_log.fillna(0, inplace=True)

        plt.figure(figsize=(18,10))
        sns.heatmap(df_log,
                    cmap="Spectral",
                    cbar_kws={"label": "-log10(p-valor)"},
                    vmin=0, vmax=30)

        plt.title(title)
        plt.xlabel("Fenotipos")
        plt.ylabel("Genotipos")
        plt.tight_layout()

        plt.savefig(heatmap_path, dpi=300)
        plt.close()

        print(f"Heatmap saved in: {heatmap_path}")

    except Exception as e:
        print(f"Error in generating heatmap form {matrix_path}: {e}")

for matrix_path in matrix_files: 
    base_name = os.path.basename(matrix_path).replace(".tsv", "").replace("_matrix","")
    heatmap_path = os.path.join(output_dir, base_name + "_heatmap.png")

    title = base_name.replace("_"," ")
    generate_heatmap(matrix_path, heatmap_path, title=title)

print("Process completed.")

#
#--------------------------------------------------------------------------
#
