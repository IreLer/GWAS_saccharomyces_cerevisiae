#--------------------------------------------------------------
#
# ILC Jun 2025
# TFM - Universidad Internacional de Valencia - Máster en Bioinformática
#
#--------------------------------------------------------------
#
# HYPEROSMOTIC STRESS RESISTANCE 
#
#-------------------------------------------------------------
#
import pandas as pd
import os
import sys
from collections import Counter
#
#--------------------------------------------------------------
#
# Required input files.

betapos_file = "05_analysis_results/YPDKCL2M_BETApos.tsv"
betaneg_file = "05_analysis_results/YPDNACL1M_BETAneg.tsv"
betapos2_file = "05_analysis_results/YPDNACL1M_BETApos.tsv"
osmotic_increased = "06_result_files/hyperosmotic_stress_increased_KCL_filtered.txt"
osmotic_increased2 = "06_result_files/hyperosmotic_stress_increased_NACL_filtered.txt"
osmotic_decreased = "06_result_files/hyperosmotic_stress_decreased_NACL_filtered.txt"
output_dir = "06_result_files"
log_file = "output_std_06_osmotic.txt"

# Create directory. 

os.makedirs(output_dir, exist_ok=True)

# Redirect std and err output to log file. 

sys.stdout = open (log_file, 'w')
sys.stderr = sys.stdout

#
#--------------------------------------------------------------
#
# 1. Compare significant genes with reference list - BETApos.
#

genes_beta_pos = []
with open(betapos_file) as f: 
    header = next(f)
    for line in f: 
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_pos.append(fields[10])
        else: 
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_pos = Counter(genes_beta_pos)

genes_osmotic_increased = set ()
with open (osmotic_increased) as f: 
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_osmotic_increased.add(fields[0])

common_genes_pos = beta_counts_pos.keys() & genes_osmotic_increased
print(f"Common genes between YPDKCL2M_BETApos and hyperosmotic_stress_increased: {len(common_genes_pos)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_pos):
      print(f"{gene}\t{beta_counts_pos[gene]}")

output_file_pos = os.path.join(output_dir, "genes_hyperosmotic_KCl2M_increased_counts.tsv")
with open (output_file_pos, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_pos):
        out.write(f"{gene}\t{beta_counts_pos[gene]}\n")

print(f"\nResults saved in: {output_file_pos}")

#
#--------------------------------------------------------------
#
# 2. Compare significant genes with reference list - BETApos.
#

genes_beta_pos_2 = []
with open(betapos2_file) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_pos_2.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_pos_2 = Counter(genes_beta_pos_2)

genes_osmotic_increased_2 = set ()
with open (osmotic_increased2) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_osmotic_increased_2.add(fields[0])

common_genes_pos_2 = beta_counts_pos_2.keys() & genes_osmotic_increased_2
print(f"Common genes between YPDNACL1M_BETApos and hyperosmotic_stress_increased: {len(common_genes_pos)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_pos_2):
      print(f"{gene}\t{beta_counts_pos_2[gene]}")

output_file_pos_2 = os.path.join(output_dir, "genes_hyperosmotic_NaCl1M_increased_counts.tsv")
with open (output_file_pos_2, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_pos_2):
        out.write(f"{gene}\t{beta_counts_pos_2[gene]}\n")

print(f"\nResults saved in: {output_file_pos_2}")

#--------------------------------------------------------------
#
# 3. Compare significant genes with reference list - BETAneg.
#
#

genes_beta_neg = []
with open(betaneg_file) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) >10:
            genes_beta_neg.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_neg = Counter(genes_beta_neg)

genes_osmotic_decreased = set ()
with open (osmotic_decreased) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_osmotic_decreased.add(fields[0])

common_genes_neg = beta_counts_neg.keys() & genes_osmotic_decreased

print(f"Common genes between YPDNACL1M_BETAneg and hyperosmotic_stress_decreased: {len(common_genes_neg)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_neg):
      print(f"{gene}\t{beta_counts_neg[gene]}")

output_file_neg = os.path.join(output_dir, "genes_hyperosmotic_NaCl1M_decreased_counts.tsv")
with open (output_file_neg, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_neg):
        out.write(f"{gene}\t{beta_counts_neg[gene]}\n")

print(f"\nResults saved in: {output_file_neg}")

#
#--------------------------------------------------------------
#

