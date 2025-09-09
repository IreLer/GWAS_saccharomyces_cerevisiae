#--------------------------------------------------------------
#
# ILC Jun 2025
# TFM - Universidad Internacional de Valencia - Máster en Bioinformática
#
#--------------------------------------------------------------
#
# SYNTHESIS INHIBITORS RESISTANCE 
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

betaneg_file_anis = "05_analysis_results/YPDANISO10_BETAneg.tsv"
betaneg_file_6au = "05_analysis_results/YPD6AU_BETAneg.tsv"

chemical_decreased_anis = "06_result_files/chemical_resistance_decreased_anisomycin_filtered.txt"
chemical_decreased_6au = "06_result_files/chemical_resistance_decreased_azauracil_filtered.txt"

betapos_file_10 = "05_analysis_results/YPDANISO10_BETApos.tsv"
betapos_file_20 = "05_analysis_results/YPDANISO20_BETApos.tsv"
betapos_file_50 = "05_analysis_results/YPDANISO50_BETApos.tsv"

chemical_increased_anis = "06_result_files/chemical_resistance_increased_anisomycin_filtered.txt"

output_dir = "06_result_files"
log_file = "output_std_06_inhibitor.txt"

# Create directory. 

os.makedirs(output_dir, exist_ok=True)

# Redirect std and err output to log file. 

sys.stdout = open (log_file, 'w')
sys.stderr = sys.stdout

#
#--------------------------------------------------------------
#
# 1. Compare significant genes with reference list - BETAneg_ANISO10.
#

genes_beta_neg_anis = []
with open(betaneg_file_anis) as f: 
    header = next(f)
    for line in f: 
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_neg_anis.append(fields[10])
        else: 
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_neg_anis = Counter(genes_beta_neg_anis)

genes_anis_decreased = set ()
with open (chemical_decreased_anis) as f: 
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_anis_decreased.add(fields[0])

common_genes_neg_anis = beta_counts_neg_anis.keys() & genes_anis_decreased
print(f"Common genes between YPDANISO10_BETAneg and chemical_stress_decreased: {len(common_genes_neg_anis)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_neg_anis):
      print(f"{gene}\t{beta_counts_neg_anis[gene]}")

output_file_neg_anis = os.path.join(output_dir, "genes_inhibitor_anisomycin_decreased_counts.tsv")
with open (output_file_neg_anis, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_neg_anis):
        out.write(f"{gene}\t{beta_counts_neg_anis[gene]}\n")

print(f"\nResults saved in: {output_file_neg_anis}")

#
#--------------------------------------------------------------
#
# 2. Compare significant genes with reference list - BETAneg_6AU.
#

genes_beta_neg_6au = []
with open(betaneg_file_6au) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_neg_6au.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_neg_6au = Counter(genes_beta_neg_6au)

genes_6au_decreased = set ()
with open (chemical_decreased_6au) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_6au_decreased.add(fields[0])

common_genes_neg_6au = beta_counts_neg_6au.keys() & genes_6au_decreased
print(f"Common genes between YPD6AU_BETAneg and chemical_stress_decreased: {len(common_genes_neg_6au)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_neg_6au):
      print(f"{gene}\t{beta_counts_neg_6au[gene]}")

output_file_neg_6au = os.path.join(output_dir, "genes_inhibitor_6AU_decreased_counts.tsv")
with open (output_file_neg_6au, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_neg_6au):
        out.write(f"{gene}\t{beta_counts_neg_6au[gene]}\n")

print(f"\nResults saved in: {output_file_neg_6au}")


#
#--------------------------------------------------------------
#
# 3. Compare significant genes with reference list - BETApos_10ANISO.
#

genes_beta_pos_10 = []
with open(betapos_file_10) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_pos_10.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_pos_10 = Counter(genes_beta_pos_10)

genes_10_increased = set ()
with open (chemical_increased_anis) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_10_increased.add(fields[0])

common_genes_pos_10 = beta_counts_pos_10.keys() & genes_10_increased
print(f"Common genes between YPDANISO10_BETApos and chemical_stress_increased: {len(common_genes_pos_10)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_pos_10):
      print(f"{gene}\t{beta_counts_pos_10[gene]}")

output_file_pos_10 = os.path.join(output_dir, "genes_inhibitor_anisomycin10_increased_counts.tsv")
with open (output_file_pos_10, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_pos_10):
        out.write(f"{gene}\t{beta_counts_pos_10[gene]}\n")

print(f"\nResults saved in: {output_file_pos_10}")

#
#--------------------------------------------------------------
#
# 4. Compare significant genes with reference list - BETApos_20ANISO.
#

genes_beta_pos_20 = []
with open(betapos_file_20) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_pos_20.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_pos_20 = Counter(genes_beta_pos_20)

genes_20_increased = set ()
with open (chemical_increased_anis) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_20_increased.add(fields[0])

common_genes_pos_20 = beta_counts_pos_20.keys() & genes_20_increased
print(f"Common genes between YPDANISO20_BETApos and chemical_stress_increased: {len(common_genes_pos_20)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_pos_20):
      print(f"{gene}\t{beta_counts_pos_20[gene]}")

output_file_pos_20 = os.path.join(output_dir, "genes_inhibitor_anisomycin20_increased_counts.tsv")
with open (output_file_pos_20, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_pos_20):
        out.write(f"{gene}\t{beta_counts_pos_20[gene]}\n")

print(f"\nResults saved in: {output_file_pos_20}")

#
#--------------------------------------------------------------
#
# 5. Compare significant genes with reference list - BETApos_50ANISO.
#

genes_beta_pos_50 = []
with open(betapos_file_50) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_pos_50.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_pos_50 = Counter(genes_beta_pos_50)

genes_50_increased = set ()
with open (chemical_increased_anis) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_50_increased.add(fields[0])

common_genes_pos_50 = beta_counts_pos_50.keys() & genes_50_increased
print(f"Common genes between YPDANISO50_BETApos and chemical_stress_increased: {len(common_genes_pos_50)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_pos_50):
      print(f"{gene}\t{beta_counts_pos_50[gene]}")

output_file_pos_50 = os.path.join(output_dir, "genes_inhibitor_anisomycin50_increased_counts.tsv")
with open (output_file_pos_50, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_pos_50):
        out.write(f"{gene}\t{beta_counts_pos_50[gene]}\n")

print(f"\nResults saved in: {output_file_pos_50}")

#
#
