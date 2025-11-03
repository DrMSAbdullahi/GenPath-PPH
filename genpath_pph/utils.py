"""
Utility functions for the GenPath-PPH framework.

This module provides workflow functions that orchestrate operations 
using the core classes PathwayDataProcessor and GenPathHomology.

Functions:
    - extract_adjacency_matrices: Generates adjacency matrices for a list 
      of KEGG pathways using PathwayDataProcessor.
    - extract_pathway_expressions: Extracts gene expression subsets 
      corresponding to KEGG pathways and filters them to match adjacency matrices.
    - compute_betti_numbers: Computes persistent path 
      homology (PPH) Betti numbers for control and disease samples, and 
      optionally saves the results to disk.
"""

# ==========================
# Imports
# ==========================
import os                  # Handling directory paths and creating folders
import pickle              # Saving/loading Python objects (e.g., dictionaries, lists)
import numpy as np         # Numerical computations and array manipulations
from core import PathwayDataProcessor, GenPathHomology
                        # Core classes for pathway processing and persistent path homology computation


def extract_adjacency_matrices(select_path_ids, kegg_map, df_all, priority_list, kegg_files_dir=''):
    """
    Extracts adjacency matrices for selected KEGG pathways.
    
    Parameters:
        select_path_ids: list of KEGG pathway IDs to process
        kegg_map: DataFrame mapping pathways to gene symbols
        df_all: gene expression DataFrame
        priority_list: conversion priority dictionary
        kegg_files_dir: optional folder containing KEGG KGML files (default=current directory)
    
    Returns:
        Dictionary of adjacency matrices keyed by pathway ID
    """
    adj_matrices = {}
    for item in select_path_ids:
        hsa_code = item
        kgml_file = os.path.join(kegg_files_dir, f"{hsa_code}.xml") if kegg_files_dir else f"{hsa_code}.xml"
        
        data_gene_symbols = df_all.index.tolist()
        pathway_gene_symbols = kegg_map[kegg_map["PathwayID"] == hsa_code].iloc[:, 2].tolist()
        pathway_df_symbols = df_all[df_all.index.isin(pathway_gene_symbols)].index.tolist()
        
        processor = PathwayDataProcessor(kgml_file, data_gene_symbols, pathway_gene_symbols, pathway_df_symbols, priority_list)
        new_adjacency_matrix_df = processor.get_adjacency_matrix()
        adj_matrices[item] = new_adjacency_matrix_df

    return adj_matrices

def extract_pathway_expressions(select_path_ids, kegg_map, df_all, priority_list, kegg_files_dir=''):
    """
    Extracts gene expression data for selected KEGG pathways, filtering to adjacency matrix genes.
    
    Parameters:
        select_path_ids: list of KEGG pathway IDs to process
        kegg_map: DataFrame mapping pathways to gene symbols
        df_all: gene expression DataFrame
        priority_list: conversion priority dictionary
        kegg_files_dir: optional folder containing KEGG KGML files (default=current directory)
    
    Returns:
        Dictionary of gene expression DataFrames keyed by pathway ID
    """
    pathways_expressions = {}
    for item in select_path_ids:
        hsa_code = item
        kgml_file = os.path.join(kegg_files_dir, f"{hsa_code}.xml") if kegg_files_dir else f"{hsa_code}.xml"
        
        pathway_gene_symbols = kegg_map[kegg_map["PathwayID"] == hsa_code].iloc[:, 2].tolist()
        pathway_gene_expression = df_all[df_all.index.isin(pathway_gene_symbols)]
        
        processor = PathwayDataProcessor(
            kgml_file, 
            df_all.index.tolist(), 
            pathway_gene_symbols, 
            pathway_gene_expression.index.tolist(), 
            priority_list
        )
        
        new_adjacency_matrix_df = processor.get_adjacency_matrix()
        missing_symbols = pathway_gene_expression.index.difference(new_adjacency_matrix_df.index)
        
        if not missing_symbols.empty:
            pathway_gene_expression = pathway_gene_expression.drop(missing_symbols)
        
        pathways_expressions[item] = pathway_gene_expression
    
    return pathways_expressions
   
def compute_betti_numbers(select_path_ids, adj_matrices, pathways_expressions,
                          target_dimension=1, filtration_scale=1, step_size=0.01,
                          class_size=17, distance_type='1-abs-correlation',
                          save_betti_dir=None, save_edges_dir=None):
    """
    Computes persistent path homology (PPH) Betti numbers for control and disease groups for one or more KEGG pathways.

    Parameters:
        select_path_ids: list of KEGG pathway IDs
        adj_matrices: dict of adjacency matrices for pathways
        pathways_expressions: dict of expression DataFrames for pathways
        target_dimension: homology dimension to compute (default=1)
        filtration_scale: range (maximum) value of filtration (default=1)
        step_size: step size for filtration values (default=0.01)
        class_size: number of samples per group (default=17)
        distance_type: type of distance for PPH (default='1-abs-correlation')
        save_betti_dir: optional folder to save Betti numbers
        save_edges_dir: optional folder to save edges at each filtration

    Returns:
        dict: {'control': {pathway_id: betti_array}, 'disease': {pathway_id: betti_array}}
    """

    if save_betti_dir and not os.path.exists(save_betti_dir):
        os.makedirs(save_betti_dir)
    if save_edges_dir and not os.path.exists(save_edges_dir):
        os.makedirs(save_edges_dir)

    betti_results = {'control': {}, 'disease': {}}
    
    for idx, item in enumerate(select_path_ids, 1):
        print(f"Processing pathway {idx}/{len(select_path_ids)}: {item}...", flush=True)
        pathway = pathways_expressions[item]
        
        # Split control and disease
        control = np.array(pathway.iloc[:, :class_size])
        disease = np.array(pathway.iloc[:, class_size:])
        
        # Adjacency edges
        rows, cols = np.nonzero(adj_matrices[item])
        edges_adj = np.array(list(zip(rows, cols)))
        
        filtration = np.arange(0, filtration_scale, step_size)
        
        # Compute Betti numbers
        betti_c, edges_c = GenPathHomology().persistent_path_homology_from_digraph(
            control, edges_adj, target_dimension, filtration=filtration, distance_type=distance_type
        )
        betti_d, edges_d = GenPathHomology().persistent_path_homology_from_digraph(
            disease, edges_adj, target_dimension, filtration=filtration, distance_type=distance_type
        )
        
        betti_results['control'][item] = betti_c
        betti_results['disease'][item] = betti_d
        
        # Saving if directories specified or default to cwd
        if save_betti_dir or save_edges_dir:
            save_dir_betti = save_betti_dir or os.getcwd()
            save_dir_edges = save_edges_dir or os.getcwd()
            
            np.save(os.path.join(save_dir_betti, f'dim_{target_dimension}_betti_numbers_{item}_control.npy'), betti_c, allow_pickle=True)
            np.save(os.path.join(save_dir_betti, f'dim_{target_dimension}_betti_numbers_{item}_disease.npy'), betti_d, allow_pickle=True)
            
            with open(os.path.join(save_dir_edges, f'dim_{target_dimension}_edges_at_all_filtration_{item}_control.pkl'), 'wb') as f:
                pickle.dump(edges_c, f)
            with open(os.path.join(save_dir_edges, f'dim_{target_dimension}_edges_at_all_filtration_{item}_disease.pkl'), 'wb') as f:
                pickle.dump(edges_d, f)

    return betti_results

