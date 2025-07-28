#--------------------------------------------------------------
#
# ILC Jun 2025
# TFM - Universidad Internacional de Valencia - Máster en Bioinformática
#
#--------------------------------------------------------------
#
# VARIANT ANNOTATION
# 35 files with significant variants.
# Reference genome: S288C.
#
#-------------------------------------------------------------
#
import gffutils
import pandas as pd
import os
import glob
import sys
#
#--------------------------------------------------------------
#
# Required input files.

gff_file = ("./ncbi_dataset/GCF_000146045.2_R64_genomic.gff")
variant_files = glob.glob("./03_filtering_results/result_ph_*_filtered.assoc.linear") 
output_dir = "04_annotation_results"
log_file = "output_std_04.txt"

# Create directory. 

os.makedirs(output_dir, exist_ok=True)

# Redirect std and err output to log file.

sys.stdout = open (log_file, 'w')
sys.stderr = sys.stdout

#
#--------------------------------------------------------------
#
# 1. Mapping chromosome nimbers to GFF chromosome names. 

chr_map = {}

with open(gff_file) as f:
    for line in f: 
        if line.startswith("##sequence-region"):
            parts = line.strip().split()
            gff_id = parts[1]
            chr_num = str(len(chr_map)+1)
            chr_map[chr_num]=gff_id

## print(map_chr)

#
#--------------------------------------------------------------
#
# 2. Generate GFF database from reference genome if not exist.

db_file = "yeast_gff.db"

if not os.path.exists(db_file):     # Create DB if not exist.
    print("Creating GFF database with gffutils...")
    gffutils.create_db(gff_file, dbfn = db_file, force=True,
                       keep_order=True, merge_strategy='merge', 
                       sort_attribute_values=True)

else:
    print("GFF database already exists.")

print("Loading database...")

db = gffutils.FeatureDB(db_file, keep_order=True)

#
#--------------------------------------------------------------
#
# 3. Find the closest gene to a variant based on its position in the genome. 

    # chr_num: chromosome number.
    # pos: variant position on the chromosome.
    # db: preloaded genomic database (gffutils). 
    # chr_map: dictionary mapping chromosome number to GFF chromosome names. 
    # Returns: a pd.Series with gene ID, gene name and product (if annotated)

def find_closest_gene(chr_num, pos, db, chr_map):
    chrom = chr_map.get(str(chr_num))       # Get GFF chr name.
    if chrom is None: 
        return pd.Series(['NA','NA','NA'])
    try:                                    # Set 100bp along the gene.
        start = max(pos -100,1)
        end = pos + 100
                                                          # Search for genes into defined region.
        region = db.region(region=(chrom, start, end),completely_within=False)
        genes = [feature for feature in region if feature.featuretype == 'gene']     # Only genes characteristics. 
        
        if genes: 
            gene = genes[0]                                              # First gen found.
            gene_id = gene.attributes.get('gene', [''])[0]               # Extract gen atributes.
            product = gene.attributes.get('product', ['NA'])[0] 
            return pd.Series([gene.id, gene_id, product])
        else: 
            return pd.Series(['NA','NA','NA'])
    except Exception as e: 
        print(f"Error processing {chrom}:{pos} -- {e}")
        return pd.Series(['NA','NA','NA'])

#   
#--------------------------------------------------------------
#
# 4. Load filtered variant files and annotate with gene information. 

for file in variant_files:

    print(f"Annotating file: {file}")
    df = pd.read_csv(file, sep="\t")
    df[['gff_id','gene_name','product']]=df.apply (lambda row: find_closest_gene (row['CHR'], row['BP'], db, chr_map), axis=1)
    
    file_name = os.path.basename(file).replace("_filtered.assoc.linear","_annotated.tsv")
    output_file = os.path.join(output_dir, file_name)

    df.to_csv(output_file, sep="\t", index=False)
    print (f"Annotated file saved: {output_file}.")

#
#--------------------------------------------------------------
#

