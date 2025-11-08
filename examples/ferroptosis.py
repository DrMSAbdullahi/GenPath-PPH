"""
Compute persistent path homology (PPH) Betti numbers for the ferroptosis pathway (hsa04216)
using real RNA-seq expression data and KEGG pathway gene mapping.

Steps:
1. Load and preprocess expression data and KEGG mappings.
2. Extract pathway-specific gene expression matrix.
3. Obtain adjacency matrix from KEGG KGML.
4. Compute Betti numbers using PPH.

Required files:
- PBMC_RNASeq_PRNA739257_log2TPM_symbol.csv
- Genes_Map_KEGG_pathways.txt
- kegg_pathways_info_clean_data_symbol.csv
- hsa04216.xml
- priority_list.pkl
"""

import pandas as pd
import numpy as np
import pickle
from core import PathwayDataProcessor, GenPathHomology
from utils import extract_pathway_expressions, extract_adjacency_matrices


# -----------------------------
# Helper functions
# -----------------------------
def get_pathway_genes(kegg_map, pathway_id):
    """Return list of genes mapped to a given KEGG pathway ID."""
    return kegg_map.loc[kegg_map["PathwayID"] == pathway_id, "Symbol"].tolist()


def extract_pathway_data(data, kegg_map, pathway_id):
    """Return subset of expression data containing only genes in the given pathway."""
    genes = get_pathway_genes(kegg_map, pathway_id)
    return data[data.index.isin(genes)]


# -----------------------------
# 1. Define inputs
# -----------------------------
hsacode = "hsa04216"
kgml_file = f"{hsacode}.xml"
path_id = [hsacode]

# Expression data
df_tpm_rna_log2 = pd.read_csv("PBMC_RNASeq_PRNA739257_log2TPM_symbol.csv")
df_tpm_rna_log2 = df_tpm_rna_log2.rename(columns={"Unnamed: 0": "Symbol"})
df_tpm_rna_log2 = df_tpm_rna_log2.set_index("Symbol")
expression_df = df_tpm_rna_log2

# KEGG mapping and filtered info
kegg_map = pd.read_csv("Genes_Map_KEGG_pathways.txt", delimiter="\t")
kegg_filtered = pd.read_csv("kegg_pathways_info_clean_data_symbol.csv")

# Priority list (open from .pkl)
with open("priority_list.pkl", "rb") as f:
    priority_list = pickle.load(f)

# Gene symbols
data_gene_symbols = expression_df.index.tolist()
pathway_expressions = extract_pathway_expressions(path_id, kegg_map, expression_df, priority_list)

pathway_gene_symbols = get_pathway_genes(kegg_map, hsacode)
pathway_df_symbols = pathway_expressions[hsacode].index.tolist()

# -----------------------------
# 2. Process pathway adjacency
# -----------------------------
processor = PathwayDataProcessor(
    kgml_file=kgml_file,
    data_gene_symbols=data_gene_symbols,
    pathway_gene_symbols=pathway_gene_symbols,
    pathway_df_symbols=pathway_df_symbols,
    priority_list=priority_list
)
adj_matrix = processor.get_adjacency_matrix()

# Alternative (direct utility)
adj_matrices = extract_adjacency_matrices(path_id, kegg_map, expression_df, priority_list)

# -----------------------------
# 3. Compute Betti numbers
# -----------------------------
pph = GenPathHomology()
betti_results_0 = pph.compute_betti_numbers(
    select_path_ids=path_id,
    adj_matrices=adj_matrices,
    pathways_expressions=pathway_expressions,
    target_dimension=0,
    filtration_scale=1,
    step_size=0.01,
    class_size=17,
    distance_type="1-abs-correlation",
    save_betti_dir=None,
    save_edges_dir=None,
)

betti_results_1 = pph.compute_betti_numbers(
    select_path_ids=path_id,
    adj_matrices=adj_matrices,
    pathways_expressions=pathway_expressions,
    target_dimension=1,
    filtration_scale=1,
    step_size=0.01,
    class_size=17,
    distance_type="1-abs-correlation",
    save_betti_dir=None,
    save_edges_dir=None,
)

# -----------------------------
# 4. Output results 
# -----------------------------
for dim in [0, 1]:
    # Select the correct results variable
    betti_results = betti_results_0 if dim == 0 else betti_results_1

    print(f"\n=== Betti Numbers (Dimension {dim}) ===")
    for cls in ["control", "disease"]:
        if cls in betti_results and hsacode in betti_results[cls]:
            betti_values = betti_results[cls][hsacode]
            # Convert np.float64 → float
            betti_values = [
                int(x) if float(x).is_integer() else round(float(x), 4)
                for x in np.ravel(betti_values)
            ]
            print(f"{cls.capitalize()} - dim {dim}: {betti_values}")
        else:
            print(f"{cls.capitalize()} - dim {dim}: Not available")

