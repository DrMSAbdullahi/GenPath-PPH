"""
run_analysis.py

Main script to compute persistent path homology (PPH) Betti numbers
for multiple KEGG pathways using RNA-seq expression data.

Author: Muhammad Sirajo Abdullahi
Date: 2025-02-14
"""

# ==========================
# Imports
# ==========================
import os                  # Handling directory paths and creating folders
import pickle              # Saving/loading Python objects (e.g., dictionaries, lists)
import pandas as pd        # DataFrame operations for expression and KEGG data
import numpy as np         # Numerical computations and array manipulations
from tqdm import tqdm      # Progress bar for loops
from utils import extract_adjacency_matrices, extract_pathway_expressions, compute_betti_numbers
                          # Utility functions for GenPath-PPH workflow

# ==========================
# Load RNA-seq expression data
# ==========================
print("Loading and preparing expression data...", flush=True)
df_tpm_rna_log2 = pd.read_csv('PBMC_RNASeq_PRNA739257_log2TPM_symbol.csv', header=0)
df_tpm_rna_log2 = df_tpm_rna_log2.rename(columns={'Unnamed: 0': 'Symbol'})
df_tpm_rna_log2 = df_tpm_rna_log2.set_index('Symbol')
df_all = df_tpm_rna_log2

# ==========================
# Load KEGG mapping data
# ==========================
print("Loading and preparing KEGG mapping data...", flush=True)
kegg_map = pd.read_csv('Genes_Map_KEGG_pathways.txt', delimiter='\t')

# ==========================
# Load selected KEGG pathway IDs
# ==========================
print("Loading selected KEGG pathways...", flush=True)
select_ids = [
    'hsa00400', 'hsa00290', 'hsa00514', 'hsa00511', 'hsa00533', 'hsa00524', 'hsa03010', 'hsa03008', 'hsa03040',
    'hsa03430', 'hsa03410', 'hsa03420', 'hsa03450', 'hsa03320', 'hsa03266', 'hsa03267', 'hsa03060', 'hsa03082',
    'hsa03083', 'hsa03260', 'hsa03265', 'hsa03264', 'hsa03020', 'hsa03013', 'hsa03022', 'hsa03050', 'hsa02010',
    'hsa01040', 'hsa03030', 'hsa04142', 'hsa04120', 'hsa04130', 'hsa04814', 'hsa04966', 'hsa04977', 'hsa04974',
    'hsa04640', 'hsa04973', 'hsa04146', 'hsa04260', 'hsa04976', 'hsa04216', 'hsa04972', 'hsa03460', 'hsa04742',
    'hsa04978', 'hsa04964', 'hsa00563', 'hsa04970', 'hsa00970', 'hsa00780', 'hsa04672', 'hsa03250', 'hsa00910',
    'hsa04064', 'hsa04060', 'hsa00470', 'hsa04144', 'hsa04145', 'hsa03018', 'hsa04610', 'hsa04080', 'hsa04141',
    'hsa03440', 'hsa04975', 'hsa04110', 'hsa04721', 'hsa04613', 'hsa04962', 'hsa04630', 'hsa00232', 'hsa04979',
    'hsa04061', 'hsa00534', 'hsa04913', 'hsa00513', 'hsa04623', 'hsa04911', 'hsa04925', 'hsa04136', 'hsa03015',
    'hsa04744', 'hsa04714', 'hsa04115', 'hsa04960', 'hsa04392', 'hsa00510', 'hsa04924', 'hsa00592', 'hsa00190',
    'hsa00630', 'hsa00440', 'hsa04122', 'hsa04961', 'hsa04614', 'hsa04066', 'hsa00130', 'hsa00040', 'hsa04622',
    'hsa04750', 'hsa04971', 'hsa04659', 'hsa04621', 'hsa04668', 'hsa04217', 'hsa00785', 'hsa00531', 'hsa00360',
    'hsa00340', 'hsa00120', 'hsa00603', 'hsa04114', 'hsa04915', 'hsa04914', 'hsa00920', 'hsa00053', 'hsa04215',
    'hsa04150', 'hsa00532', 'hsa04918', 'hsa00515', 'hsa00730', 'hsa04919', 'hsa04330', 'hsa04140', 'hsa04512',
    'hsa00860', 'hsa00062', 'hsa00750', 'hsa04929', 'hsa04662', 'hsa04625', 'hsa00220', 'hsa04350', 'hsa04520',
    'hsa04922', 'hsa04611', 'hsa04218', 'hsa04012', 'hsa04137', 'hsa04152', 'hsa04923', 'hsa04540', 'hsa04514',
    'hsa04927', 'hsa00100', 'hsa04620', 'hsa00052', 'hsa04928', 'hsa04380', 'hsa04270', 'hsa04390', 'hsa00061',
    'hsa04666', 'hsa00740', 'hsa04022', 'hsa04710', 'hsa04071', 'hsa04360', 'hsa00770', 'hsa00450', 'hsa04211',
    'hsa04340', 'hsa04920', 'hsa04664', 'hsa00260', 'hsa04912', 'hsa00650', 'hsa00790', 'hsa04650', 'hsa04370',
    'hsa04658', 'hsa04072', 'hsa04921', 'hsa04210', 'hsa04910', 'hsa04935', 'hsa00900', 'hsa04917', 'hsa04726',
    'hsa00604', 'hsa00983', 'hsa04660', 'hsa00640', 'hsa04722', 'hsa00380', 'hsa00250', 'hsa00330', 'hsa04550',
    'hsa04213', 'hsa00270', 'hsa00020', 'hsa04024', 'hsa04916', 'hsa00310', 'hsa04510', 'hsa04728', 'hsa00520',
    'hsa00350', 'hsa04724', 'hsa04730', 'hsa04068', 'hsa04530', 'hsa00010', 'hsa04810', 'hsa00280', 'hsa00430',
    'hsa04657', 'hsa04020', 'hsa04670', 'hsa04151', 'hsa04725', 'hsa04720', 'hsa00620', 'hsa00410', 'hsa00071',
    'hsa04926', 'hsa00051', 'hsa00982', 'hsa00512', 'hsa04015', 'hsa04723', 'hsa04612', 'hsa04727', 'hsa04261',
    'hsa00601', 'hsa00500', 'hsa00030', 'hsa04310', 'hsa04014', 'hsa00591', 'hsa04010', 'hsa04371', 'hsa04713',
    'hsa04740', 'hsa00140', 'hsa00830', 'hsa00670', 'hsa04062', 'hsa00590', 'hsa00561', 'hsa00480', 'hsa00760',
    'hsa00565', 'hsa00562', 'hsa00240', 'hsa00564', 'hsa00600', 'hsa00980', 'hsa04070', 'hsa00230'
]

# ==========================
# Load gene conversion priority list
# ==========================
print("Loading gene conversion priority list...", flush=True)
with open('kegg_genes_ids_to_symbols_priority_list.pkl', 'rb') as f:
    priority_list = pickle.load(f)

# ==========================
# Extract adjacency matrices and pathway expressions
# ==========================
print("Extracting adjacency matrices and pathway expressions...", flush=True)
adj_matrices = extract_adjacency_matrices(select_ids, kegg_map, df_all, priority_list)
pathways_expressions = extract_pathway_expressions(select_ids, kegg_map, df_all, priority_list)

# ==========================
# Check consistency of dimensions
# ==========================
print("Checking adjacency matrix vs. expression dimensions...", flush=True)
adj_length = [adj.shape[0] for adj in adj_matrices.values()]
expression_length = [exp.shape[0] for exp in pathways_expressions.values()]
if adj_length != expression_length:
    raise ValueError("Mismatch: adjacency matrix and expression data lengths do not match!")

# ==========================
# Choose target dimensions to compute
# ==========================
target_dimensions = [0, 1]  # Example: compute 0-dim and 1-dim Betti numbers

# ==========================
# Compute Betti numbers for all pathways and dimensions
# ==========================
for dim in target_dimensions:
    print(f"\nComputing Betti numbers for target dimension {dim}...", flush=True)
    
    # Use tqdm to show progress per pathway
    betti_results_dim = {}
    for pathway_id in tqdm(select_ids, desc=f"Dimension {dim}"):
        result = compute_betti_numbers(
            select_path_ids=[pathway_id],
            adj_matrices=adj_matrices,
            pathways_expressions=pathways_expressions,
            target_dimension=dim,
            filtration_scale=1,        # 1 for '1-abs-correlation', 2 for '1-correlation'
            step_size=0.01,            # default for GenPath-PPH
            save_betti_dir=None,       # saves in current directory if None
            save_edges_dir=None
        )
        betti_results_dim[pathway_id] = result

print("\nDone computing Betti numbers for all selected KEGG pathways and dimensions.", flush=True)

