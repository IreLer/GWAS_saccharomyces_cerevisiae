#--------------------------------------------------------------
#
# ILC Jun 2025
# TFM - Universidad Internacional de Valencia - Máster en Bioinformática
#
#--------------------------------------------------------------
#
# OXIDATIVE AND METAL STRESS RESISTANCE 
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

betapos_file_ox = "05_analysis_results/YPDSODIUMMETAARSENITE_BETApos.tsv"
betaneg_file_ox = "05_analysis_results/YPDSODIUMMETAARSENITE_BETAneg.tsv"
oxidative_increased = "06_result_files/oxidative_stress_resistance_increased_paraquat_filtered.txt"
oxidative_decreased = "06_result_files/oxidative_stress_resistance_decreased_paraquat_filtered.txt"

betapos_file_met = "05_analysis_results/YPDCUSO410MM_BETApos.tsv"
betaneg_file_met = "05_analysis_results/YPDCUSO410MM_BETAneg.tsv"
metal_increased = "06_result_files/metal_resistance_increased_copper_filtered.txt"
metal_decreased = "06_result_files/metal_resistance_decreased_copper_filtered.txt"

output_dir = "06_result_files"
log_file = "output_std_06_oxidative.txt"

# Create directory. 

os.makedirs(output_dir, exist_ok=True)

# Redirect std and err output to log file. 

sys.stdout = open (log_file, 'w')
sys.stderr = sys.stdout

#
#--------------------------------------------------------------
#
# 1. Compare significant genes with reference list - oxidative_BETApos.
#

genes_beta_pos_ox = []
with open(betapos_file_ox) as f: 
    header = next(f)
    for line in f: 
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_pos_ox.append(fields[10])
        else: 
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_pos_ox = Counter(genes_beta_pos_ox)

genes_oxidative_increased = set ()
with open (oxidative_increased) as f: 
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_oxidative_increased.add(fields[0])

common_genes_pos_ox = beta_counts_pos_ox.keys() & genes_oxidative_increased
print(f"Common genes between YPDSODIUMMARSENITE_BETApos and oxidative_stress_increased: {len(common_genes_pos_ox)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_pos_ox):
      print(f"{gene}\t{beta_counts_pos_ox[gene]}")

output_file_pos_ox = os.path.join(output_dir, "genes_oxidative_paraquat_increased_counts.tsv")
with open (output_file_pos_ox, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_pos_ox):
        out.write(f"{gene}\t{beta_counts_pos_ox[gene]}\n")

print(f"\nResults saved in: {output_file_pos_ox}")

#
#--------------------------------------------------------------
#
# 2. Compare significant genes with reference list - oxidative_BETAneg.
#
#

genes_beta_neg_ox = []
with open(betaneg_file_ox) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) >10:
            genes_beta_neg_ox.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_neg_ox = Counter(genes_beta_neg_ox)

genes_oxidative_decreased = set ()
with open (oxidative_decreased) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_oxidative_decreased.add(fields[0])

common_genes_neg_ox = beta_counts_neg_ox.keys() & genes_oxidative_decreased

print(f"Common genes between YPDSODIUMMETAARSENITE_BETAneg and oxidative_stress_decreased: {len(common_genes_neg_ox)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_neg_ox):
      print(f"{gene}\t{beta_counts_neg_ox[gene]}")

output_file_neg_ox = os.path.join(output_dir, "genes_oxidative_paraquat_decreased_counts.tsv")
with open (output_file_neg_ox, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_neg_ox):
        out.write(f"{gene}\t{beta_counts_neg_ox[gene]}\n")

print(f"\nResults saved in: {output_file_neg_ox}")

#
#--------------------------------------------------------------
#
# 3. Compare significant genes with reference list - metal_BETApos.
#

genes_beta_pos_met = []
with open(betapos_file_met) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 10:
            genes_beta_pos_met.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_pos_met = Counter(genes_beta_pos_met)

genes_metal_increased = set ()
with open (metal_increased) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_metal_increased.add(fields[0])

common_genes_pos_met = beta_counts_pos_met.keys() & genes_metal_increased
print(f"Common genes between YPDCUSO410MM_BETApos and metal_stress_increased: {len(common_genes_pos_met)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_pos_met):
      print(f"{gene}\t{beta_counts_pos_met[gene]}")

output_file_pos_met = os.path.join(output_dir, "genes_metal_copper_increased_counts.tsv")
with open (output_file_pos_met, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_pos_met):
        out.write(f"{gene}\t{beta_counts_pos_met[gene]}\n")

print(f"\nResults saved in: {output_file_pos_met}")

#
#--------------------------------------------------------------
#
# 4. Compare significant genes with reference list - metal_BETAneg.
#
#

genes_beta_neg_met = []
with open(betaneg_file_met) as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) >10:
            genes_beta_neg_met.append(fields[10])
        else:
            print(f"Line with less than 11 columns skipped: {fields}")

beta_counts_neg_met = Counter(genes_beta_neg_met)

genes_metal_decreased = set ()
with open (metal_decreased) as f:
    for line in f:
        if not line.startswith("!"):
            fields = line.strip().split("\t")
            genes_metal_decreased.add(fields[0])

common_genes_neg_met = beta_counts_neg_met.keys() & genes_metal_decreased

print(f"Common genes between YPDCUSO410MM_BETAneg and metal_stress_decreased: {len(common_genes_neg_met)}\n")
print(f"Gene\tCounts")
for gene in sorted(common_genes_neg_met):
      print(f"{gene}\t{beta_counts_neg_met[gene]}")

output_file_neg_met = os.path.join(output_dir, "genes_metal_copper_decreased_counts.tsv")
with open (output_file_neg_met, "w") as out:
      out.write("Gene\tCounts\n")
      for gene in sorted(common_genes_neg_met):
        out.write(f"{gene}\t{beta_counts_neg_met[gene]}\n")

print(f"\nResults saved in: {output_file_neg_met}")

#
#--------------------------------------------------------------
#
