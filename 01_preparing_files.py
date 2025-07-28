#--------------------------------------------------------------
#
# ILC Jun 2025
# TFM - Universidad Internacional de Valencia - Máster en Bioinformática
#
#--------------------------------------------------------------
#
# GWAS FILE PREPARATION
# 1011Matrix.gvcf.gz -- Variant matrix file.
# phenoMatrix_35ConditionsNormalizedByYPD.tab.gz -- Phenotype table file.
#
#--------------------------------------------------------------
#
import subprocess
import pandas as pd
import os
import sys
import gzip
#
#--------------------------------------------------------------
#
# Required input files. 

gvcf_path = "1011Matrix.gvcf.gz"
pheno_path = "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz"
log_file = "output_std_01.txt"
output_dir = "01_modified_files"
vcf_path = os.path.join(output_dir, "1011Matrix.vcf")
pheno_dir = "01_modified_files/phenotype_files"

# Create directory.

os.makedirs(output_dir, exist_ok=True)
os.makedirs(pheno_dir, exist_ok=True)

# Redirect std and err output to log file. 

sys.stdout = open (log_file, 'w') 
sys.stderr = sys.stdout

#
#-------------------------------------------------------------
#
# 1. Generate .vcf file from .gvcf
#

def index_gvcf (gvcf_path):
    print(f"Indexing {gvcf_path} with GATK...")
    cmd_1 = ["gatk", "IndexFeatureFile", "-I", gvcf_path]
    subprocess.run(cmd_1, check=True)
    print("Indexing completed.")

def gvcf_to_vcf(gvcf_path, vcf_path):
    print(f"Converting {gvcf_path} to {vcf_path} with GATK...")
    cmd_2 = [ "gatk", "SelectVariants","-V", gvcf_path,
           "-O", vcf_path, "--exclude-non-variants"]
    subprocess.run(cmd_2, check=True)
    print("Conversion completed.")

index_gvcf (gvcf_path)
gvcf_to_vcf(gvcf_path, vcf_path)

#
#--------------------------------------------------------------
#
# 2. Rename chr in .vcf file.
# 

def vcf_change (vcf_input, vcf_output):
    print(f"Fixing chromosome names...")
    cmd_3 = f"sed 's/^chromosome//' {vcf_input} > {vcf_output}"
    subprocess.run(cmd_3, shell=True, check=True)
    print("Chromosome names fixed successfully.")
    
vcf_input = "01_modified_files/1011Matrix.vcf"
vcf_output = "01_modified_files/1011Matrix_fixed_chr.vcf"

vcf_change (vcf_input, vcf_output)

#
#--------------------------------------------------------------
#
# 3. Generate .bed, .bim and .fam files.
#

def plink_make_bed (vcf_entry, prefix_output):
    print(f"Generating required files...")
    cmd_4 = ["plink", "--vcf", vcf_entry, "--biallelic-only", "strict", 
             "--maf", "0.01", "--make-bed", "--out", prefix_output]
    subprocess.run(cmd_4, check=True)
    print("Required files generated successfully.")

vcf_entry = "01_modified_files/1011Matrix_fixed_chr.vcf"
prefix_output = "01_modified_files/GWAS_Matrix"

plink_make_bed (vcf_entry, prefix_output)

#
#--------------------------------------------------------------
#
# 4. Phenotype file modification. 
#
    # Remove 'SACE_' prefix from pheno file.

def clean_prefix (input_path, clean_path):
    print(f"Removing 'SACE_' prefix from {input_path}...")
    with gzip.open(input_path, 'rt') as infile, open(clean_path, 'w') as temp_file:
        for line in infile:
            clean_line = line.replace("SACE_","")
            temp_file.write(clean_line)
    print("All 'SACE_' occurrences removed.")

pheno_clean = "01_modified_files/pheno_clean.tab"

clean_prefix(pheno_path,pheno_clean)

    # Duplicate first column and name them as FID and IID.

pheno_FID_IID = "01_modified_files/phenoMatrix_FID_IID.tab"

df = pd.read_csv(pheno_clean, sep='\t', index_col=0)

df.insert(0,'IID', df.index)
df.insert(0,'FID', df.index)

df.to_csv(pheno_FID_IID, sep='\t', index=False)

  # Phenotype names.

df_pheno = pd.read_csv(pheno_FID_IID, sep='\t')
ph = df_pheno.columns.tolist()

    # Unique file for each ph.

for i in range(2,37):
    col_name = ph[i]
    print(f"Generating file for: {col_name}...")
    df_tmp = df_pheno.iloc[:, [0,1,i]]
    final_path = os.path.join(pheno_dir, f"pheno_{col_name}.txt")
    df_tmp.to_csv(final_path, sep='\t', header=True, index=False)
    print(f"File {final_path} successfully generated.")

print(f"Completed: all files have been succesfully generated.")

#
#--------------------------------------------------------------
#
