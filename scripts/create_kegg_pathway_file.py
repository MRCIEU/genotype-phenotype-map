import gzip
import os

data_dir <- Sys.getenv("DATA_DIR")
variant_annotation_dir = f"{data_dir}variant_annotation/"

# Your exact local file paths
local_uniprot_file = f"{variant_annotation_dir}/HUMAN_9606_idmapping.dat.gz"
local_kegg_pathway = f"{variant_annotation_dir}/kegg_pathway_map.txt"
local_pathway_descriptions = f"{variant_annotation_dir}/kegg_pathway_descriptions.txt"

output_database_file = f"{variant_annotation_dir}/kegg_pathways_output.txt"

# 1. Map UniProt Accession to GeneID (NCBI) and Ensembl (ENSG)
uniprot_to_ncbi = {}
uniprot_to_ensg = {}

print("Parsing UniProt idmapping to link Ensembl and NCBI Gene IDs via protein bridges...")

with gzip.open(local_uniprot_file, 'rt', encoding='utf-8') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
        uniprot_acc, db_type, db_value = parts[0], parts[1], parts[2]
        
        # Capture the Ensembl Gene ID (ENSG) and strip any version decimals
        if db_type == "Ensembl":
            uniprot_to_ensg[uniprot_acc] = db_value.split('.')[0]
            
        # Capture the NCBI Entrez/KEGG Gene ID
        elif db_type == "GeneID":
            uniprot_to_ncbi[uniprot_acc] = db_value

# 2. Bridge the two tracking systems: Build direct NCBI -> ENSG Map
ncbi_to_ensg = {}
for uniprot_acc, ensg_id in uniprot_to_ensg.items():
    ncbi_id = uniprot_to_ncbi.get(uniprot_acc)
    if ncbi_id:
        ncbi_to_ensg[ncbi_id] = ensg_id

print(f"Successfully bridged {len(ncbi_to_ensg)} NCBI Gene IDs directly to clean ENSG IDs.")

# 3. Read the Pathway Names into memory from pathway_descriptions.txt
print("Loading pathway descriptions...")
pathway_names = {}
with open(local_pathway_descriptions, 'r', encoding='utf-8') as f:
    for line in f:
        stripped = line.strip()
        if not stripped: 
            continue
        parts = stripped.split('\t')
        if len(parts) < 2: 
            continue
        
        # Clean up path prefixes if they exist
        clean_path = parts[0].replace("path:", "")
        clean_desc = parts[1].replace(" - Homo sapiens (human)", "")
        pathway_names[clean_path] = clean_desc

# 4. Process kegg_pathway.txt to produce the final relational layout
print("Building final table rows...")
rows_written = 0

with open(local_kegg_pathway, 'r', encoding='utf-8') as f_in, open(output_database_file, 'w', encoding='utf-8') as f_out:
    # Set up clean header
    f_out.write("gene_id\tterm_id\tsource\tdescription\n")
    
    for line in f_in:
        stripped = line.strip()
        if not stripped: 
            continue
        parts = stripped.split('\t')
        if len(parts) < 2: 
            continue
        
        term_id = parts[0].replace("path:", "")  # e.g., 'hsa00010'
        ncbi_id = parts[1].replace("hsa:", "")   # e.g., '10327'
        
        # Look up the translated ENSG using the bridged map
        ensg_id = ncbi_to_ensg.get(ncbi_id)
        description = pathway_names.get(term_id, "")
        
        if ensg_id:
            f_out.write(f"{ensg_id}\t{term_id}\tKEGG\t{description}\n")
            rows_written += 1

print(f"Done! Saved {rows_written} mapped pathway records to: {os.path.abspath(output_database_file)}")