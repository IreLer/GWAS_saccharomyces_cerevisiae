#--------------------------------------------------------------
#
# ILC Jun 2025
# TFM - Universidad Internacional de Valencia - Máster en Bioinformática
#
#--------------------------------------------------------------
#
# GENOTYPE-PHENOTYPE ASSOCIATION ANALYSIS
# 906 variants from GWAS_Matrix files.
# 35 conditions (phenotypes) from phenoMatrix file.
#
#-------------------------------------------------------------
#
import os
import sys
import subprocess
import glob
#
#--------------------------------------------------------------
#
# Required input files.

variant_file = "01_modified_files/GWAS_Matrix"
pheno_file = "01_modified_files/phenotype_files/"
results_dir = "02_GWAS_results/"
pheno_dir = "02_GWAS_results/pheno_for_GWAS/"
log_file = "output_std_02.txt"

# Create directory. 

os.makedirs(results_dir, exist_ok=True)
os.makedirs(pheno_dir, exist_ok=True)

# Redirect std and err output to log file. 

sys.stdout = open (log_file, 'w')
sys.stderr = sys.stdout

#
#--------------------------------------------------------------
#
# 1. Genotype - phenotype association analysis for each phenotype file. 

raw_files = sorted(glob.glob(os.path.join(pheno_file, "pheno_*.txt")))
for file in raw_files:
    file_name = os.path.basename(file)
    output_file = os.path.join(pheno_dir, file_name)

    with open(file, 'r') as fin, open(output_file, 'w') as fout:
        next(fin)
        for line in fin: 
            fout.write(line)

print(f"Headers removed from phenotype files")

pheno_files = sorted(glob.glob(os.path.join(pheno_dir, "pheno_*.txt")))
total_ph = len(pheno_files)
count = 0

for ph_file in pheno_files:
    
    base = os.path.basename(ph_file)
    pheno_name = base.replace ("pheno_", "").replace (".txt", "")
    count += 1

    print(f"[{count}/{total_ph}] Analizando fenotipo: {pheno_name}...")
    
    out_path = os.path.join (results_dir, f"result_ph_{pheno_name}")
    cmd = [ "plink", "--bfile", variant_file, 
            "--pheno", ph_file, "--allow-no-sex", "--linear", 
            "--out", out_path]
    
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    
    for line in process.stdout:
        print(line, end='')
    process.wait()

    if process.returncode != 0:
        print(f"error running PLINK for {pheno_name}.")

print(f"Analysis successfully completef for {total_ph} phenotypes.")

#
#-------------------------------------------------------------
#
