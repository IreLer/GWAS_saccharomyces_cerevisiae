#--------------------------------------------------------------
#
# ILC Jun 2025
# TFM - Universidad Internacional de Valencia - Máster en Bioinformática
#
#--------------------------------------------------------------
#
# RESISTANCE TO CHEMICALS 
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

betaneg_file_ben = "05_analysis_results/YPDBENOMYL500_BETAneg.tsv"
betaneg_file_caf40 = "05_analysis_results/YPDCAFEIN40_BETAneg.tsv"
betaneg_file_caf50 = "05_analysis_results/YPDCAFEIN50_BETAneg.tsv"

chemical_decreased_ben = "06_result_files/chemical_resistance_decreased_benomyl_filtered.txt"
chemical_decreased_caf = "06_result_files/chemical_resistance_decreased_caffeine_filtered.txt"

betapos_file_ben = "05_analysis_results/YPDBENOMYL500_BETApos.tsv"
betapos_file_caf40 = "05_analysis_results/YPDCAFEIN40_BETApos.tsv"
betapos_file_caf50 = "05_analysis_results/YPDCAFEIN50_BETApos.tsv"

chemical_increased_ben = "06_result_files/chemical_resistance_increased_benomyl_filtered.txt"
chemical_increased_caf = "06_result_files/chemical_resistance_increased_caffeine_filtered.txt"

output_dir = "06_result_files"
log_file = "output_std_06_chemical.txt"

# Create directory. 

os.makedirs(output_dir, exist_ok=True)

# Redirect std and err output to log file. 

sys.stdout = open (log_file, 'w')
sys.stderr = sys.stdout

#
#--------------------------------------------------------------
#
# 1. Compare significant genes with reference list - BETAneg_BENOMYL.
#

genes_beta_neg_ben = []
with open(betaneg_file_ben) as f: 
    header = next(f)
    for line in f: 
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_neg_ben.append(fields[10])
        else: 
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_neg_ben = Counter(genes_beta_neg_ben)

genes_ben_decreased = set ()
with open (chemical_decreased_ben) as f: 
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_ben_decreased.add(fields[0])

common_genes_neg_ben = beta_counts_neg_ben.keys() & genes_ben_decreased
print(f"Common genes between YPDBENOMYL200_BETAneg and chemical_stress_decreased: {len(common_genes_neg_ben)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_neg_ben):
      print(f"{gene}\t{beta_counts_neg_ben[gene]}")

output_file_neg_ben = os.path.join(output_dir, "genes_chemical_benomyl_decreased_counts.tsv")
with open (output_file_neg_ben, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_neg_ben):
        out.write(f"{gene}\t{beta_counts_neg_ben[gene]}\n")

print(f"\nResults saved in: {output_file_neg_ben}")

#
#--------------------------------------------------------------
#
# 2. Compare significant genes with reference list - BETAneg_CAFFEINE40.
#

genes_beta_neg_caf40 = []
with open(betaneg_file_caf40) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_neg_caf40.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_neg_caf40 = Counter(genes_beta_neg_caf40)

genes_caf40_decreased = set ()
with open (chemical_decreased_caf) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_caf40_decreased.add(fields[0])

common_genes_neg_caf40 = beta_counts_neg_caf40.keys() & genes_caf40_decreased
print(f"Common genes between YPDCAFEIN40_BETAneg and chemical_stress_decreased: {len(common_genes_neg_caf40)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_neg_caf40):
      print(f"{gene}\t{beta_counts_neg_caf40[gene]}")

output_file_neg_caf40 = os.path.join(output_dir, "genes_chemical_caffeine40_decreased_counts.tsv")
with open (output_file_neg_caf40, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_neg_caf40):
        out.write(f"{gene}\t{beta_counts_neg_caf40[gene]}\n")

print(f"\nResults saved in: {output_file_neg_caf40}")

#
#--------------------------------------------------------------
#
# 3. Compare significant genes with reference list - BETAneg_CAFFEINE50.
#

genes_beta_neg_caf50 = []
with open(betaneg_file_caf50) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_neg_caf50.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_neg_caf50 = Counter(genes_beta_neg_caf50)

genes_caf50_decreased = set ()
with open (chemical_decreased_caf) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_caf50_decreased.add(fields[0])

common_genes_neg_caf50 = beta_counts_neg_caf50.keys() & genes_caf50_decreased
print(f"Common genes between YPDCAFEIN50_BETAneg and chemical_stress_decreased: {len(common_genes_neg_caf50)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_neg_caf50):
      print(f"{gene}\t{beta_counts_neg_caf50[gene]}")

output_file_neg_caf50 = os.path.join(output_dir, "genes_chemical_caffeine50_decreased_counts.tsv")
with open (output_file_neg_caf50, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_neg_caf50):
        out.write(f"{gene}\t{beta_counts_neg_caf50[gene]}\n")

print(f"\nResults saved in: {output_file_neg_caf50}")

#
#--------------------------------------------------------------
#
# 4. Compare significant genes with reference list - BETApos_BENOMYL.
#

genes_beta_pos_ben = []
with open(betapos_file_ben) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_pos_ben.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_pos_ben = Counter(genes_beta_pos_ben)

genes_ben_increased = set ()
with open (chemical_increased_ben) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_ben_increased.add(fields[0])

common_genes_pos_ben = beta_counts_pos_ben.keys() & genes_ben_increased
print(f"Common genes between YPDBENOMYL200_BETApos and chemical_stress_increased: {len(common_genes_pos_ben)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_pos_ben):
      print(f"{gene}\t{beta_counts_pos_ben[gene]}")

output_file_pos_ben = os.path.join(output_dir, "genes_chemical_benomyl_increased_counts.tsv")
with open (output_file_pos_ben, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_pos_ben):
        out.write(f"{gene}\t{beta_counts_pos_ben[gene]}\n")

print(f"\nResults saved in: {output_file_pos_ben}")

#
#--------------------------------------------------------------
#
# 5. Compare significant genes with reference list - BETApos_CAFFEINE40.
#

genes_beta_pos_caf40 = []
with open(betapos_file_caf40) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_pos_caf40.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_pos_caf40 = Counter(genes_beta_pos_caf40)

genes_caf40_increased = set ()
with open (chemical_increased_caf) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_caf40_increased.add(fields[0])

common_genes_pos_caf40 = beta_counts_pos_caf40.keys() & genes_caf40_increased
print(f"Common genes between YPDCAFEIN40_BETApos and chemical_stress_increased: {len(common_genes_pos_caf40)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_pos_caf40):
      print(f"{gene}\t{beta_counts_pos_caf40[gene]}")

output_file_pos_caf40 = os.path.join(output_dir, "genes_chemical_caffeine40_increased_counts.tsv")
with open (output_file_pos_caf40, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_pos_caf40):
        out.write(f"{gene}\t{beta_counts_pos_caf40[gene]}\n")

print(f"\nResults saved in: {output_file_pos_caf40}")

#
#--------------------------------------------------------------
#
# 6. Compare significant genes with reference list - BETApos_CAFFEINE50.
#

genes_beta_pos_caf50 = []
with open(betapos_file_caf50) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_pos_caf50.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_pos_caf50 = Counter(genes_beta_pos_caf50)

genes_caf50_increased = set ()
with open (chemical_increased_caf) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_caf50_increased.add(fields[0])

common_genes_pos_caf50 = beta_counts_pos_caf50.keys() & genes_caf50_increased
print(f"Common genes between YPDCAFEIN50_BETApos and chemical_stress_increased: {len(common_genes_pos_caf50)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_pos_caf50):
      print(f"{gene}\t{beta_counts_pos_caf50[gene]}")

output_file_pos_caf50 = os.path.join(output_dir, "genes_chemical_caffeine50_increased_counts.tsv")
with open (output_file_pos_caf50, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_pos_caf50):
        out.write(f"{gene}\t{beta_counts_pos_caf50[gene]}\n")

print(f"\nResults saved in: {output_file_pos_caf50}")

#
#--------------------------------------------------------------
#
