"""
Utility functions for the GenPath-PPH framework.

This module provides workflow functions that orchestrate operations 
using the core classes PathwayDataProcessor and GenPathHomology.

Functions:
    - extract_adjacency_matrices: Generates adjacency matrices for a list 
      of KEGG pathways using PathwayDataProcessor.
    - extract_pathway_expressions: Extracts gene expression subsets 
      corresponding to KEGG pathways and filters them to match adjacency matrices.
"""

# ==========================
# Imports
# ==========================
import os                  # Handling directory paths and creating folders
import numpy as np         # Numerical computations and array manipulations
import pandas as pd          # DataFrame operations for gene expression and adjacency matrices
import matplotlib.pyplot as plt      # Visualization of persistence landscapes and results
from core import PathwayDataProcessor, GenPathHomology
                        # Core classes for pathway processing and persistent path homology computation
from scipy.stats import ks_2samp
from statsmodels.stats.multitest import multipletests

# ------------------------------
# Data and pathway preprocessing
# ------------------------------

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

# -------------------
# Statistical testing
# -------------------

def pooled_std(control, disease):
    """
    Computes pooled standard deviation between two samples.

    Parameters:
        control: array-like, control group values
        disease: array-like, disease group values

    Returns:
        float: pooled standard deviation
    """
    n1, n2 = len(control), len(disease)
    var1, var2 = np.var(control, ddof=1), np.var(disease, ddof=1)
    return np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

def permutation_test_effect_size(betti_control, betti_disease, num_permutations=1000, effect_type="ks", seed=0):
    """
    Computes permutation-based effect size or KS statistic between two samples.

    Parameters:
        betti_control: np.array, Betti numbers for the control group
        betti_disease: np.array, Betti numbers for the disease group
        num_permutations: int, number of random permutations (default=1000)
        effect_type: str, 'cohen_d' for effect size or 'ks' for Kolmogorov-Smirnov test
        seed: int, random seed for reproducibility (default=0)

    Returns:
        tuple: (observed_effect, mean_difference, empirical_p_value)
    """
    np.random.seed(seed)  # Set the random seed for reproducibility
    
    # Compute observed effect
    if effect_type == "cohen_d":
        # Skip empty arrays
        if len(betti_control) == 0 or len(betti_disease) == 0:
            return np.nan, np.nan, np.nan

        pooled_stdev = pooled_std(betti_control, betti_disease)
        mean_diff = np.mean(betti_control) - np.mean(betti_disease)

        if pooled_stdev == 0:
            observed_effect = np.nan
        else:
            observed_effect = mean_diff / pooled_stdev

    elif effect_type == "ks":
        observed_effect = ks_2samp(betti_control, betti_disease).statistic
        mean_diff = np.nan  # not applicable for KS

    else:
        raise ValueError("Invalid effect_type. Choose 'cohen_d' or 'ks'.")

    # Combine and permute
    combined = np.concatenate([betti_control, betti_disease])
    permuted_effects = []

    for _ in range(num_permutations):
        np.random.shuffle(combined)
        perm_control = combined[:len(betti_control)]
        perm_disease = combined[len(betti_control):]

        if effect_type == "cohen_d":
            pooled_stdev_perm = pooled_std(perm_control, perm_disease)
            if pooled_stdev_perm == 0:
                perm_effect = np.nan
            else:
                perm_effect = (np.mean(perm_control) - np.mean(perm_disease)) / pooled_stdev_perm
        else:
            perm_effect = ks_2samp(perm_control, perm_disease).statistic

        permuted_effects.append(perm_effect)

    permuted_effects = np.array(permuted_effects)
    empirical_p = np.mean(np.abs(permuted_effects) >= np.abs(observed_effect))

    return observed_effect, mean_diff, empirical_p

def perform_ks_and_effectsize_tests(betti_dict, path_ids, num_permutations=5000, seed=0):
    """
    Performs permutation-based Cohen's d and Kolmogorov–Smirnov (KS) tests for multiple pathways.

    Parameters:
        betti_dict: dictionary containing Betti numbers for each pathway, with keys ['control'] and ['disease']
        path_ids: list of KEGG pathway IDs to analyze
        num_permutations: number of random permutations for significance testing (default=5000)
        seed: random seed for reproducibility (default=0)

    Returns:
        DataFrame containing:
            - path_id: pathway identifier (KEGG ID/hsacode)
            - mean_diff_observed: mean difference between disease and control Betti values
            - es_observed: observed Cohen’s d effect size
            - es_raw_pvalue: unadjusted p-value for Cohen’s d
            - ks_observed: observed KS statistic
            - ks_raw_pvalue: unadjusted p-value for KS test
            - es_corrected_pvalue: Benjamini–Hochberg adjusted p-value for Cohen’s d
            - ks_corrected_pvalue: Benjamini–Hochberg adjusted p-value for KS test
    """
    
    results = {
        'path_id': [],
        'mean_diff_observed': [],
        'es_observed': [],
        'es_raw_pvalue': [],
        'ks_observed': [],
        'ks_raw_pvalue': []
    }

    for idx, pid in enumerate(path_ids, 1):
        print(f"Testing pathway {idx}/{len(path_ids)}: {pid}")
        control = betti_dict[pid]['control']
        disease = betti_dict[pid]['disease']

        es, mean_diff, raw_p_es = permutation_test_effect_size(control, disease, num_permutations, "cohen_d", seed)
        ks_stat, _, raw_p_ks = permutation_test_effect_size(control, disease, num_permutations, "ks", seed)

        results['path_id'].append(pid)
        results['mean_diff_observed'].append(mean_diff)
        results['es_observed'].append(es)
        results['es_raw_pvalue'].append(raw_p_es)
        results['ks_observed'].append(ks_stat)
        results['ks_raw_pvalue'].append(raw_p_ks)

    # Adjust p-values (Benjamini-Hochberg)
    _, results['es_corrected_pvalue'], _, _ = multipletests(results['es_raw_pvalue'], method='fdr_bh')
    _, results['ks_corrected_pvalue'], _, _ = multipletests(results['ks_raw_pvalue'], method='fdr_bh')

    return pd.DataFrame(results)
    
def permutation_test_effect_size(betti_control, betti_disease, num_permutations=1000, effect_type="cohen_d", seed=0):
    """
    Computes the effect size between Betti number distributions of control and disease, 
    then performs a permutation test to estimate the empirical p-value.

    Parameters:
    - betti_control: np.array, Betti numbers for control group
    - betti_disease: np.array, Betti numbers for disease group
    - num_permutations: int, number of random permutations
    - effect_type: str, type of effect size ("cohen_d" or "ks")

    Returns:
    - observed_effect: Observed effect size
    - empirical_p: Proportion of permuted effect sizes more extreme than observed
    """

    np.random.seed(seed)  # Set the random seed for reproducibility
        
    # Compute observed effect size
    if effect_type == "cohen_d":
        # Skip empty arrays
        if len(betti_control) == 0 or len(betti_disease) == 0:
            observed_effect = np.nan

        pooled_stdev = pooled_std(betti_control, betti_disease)
        # Handle cases where pooled standard deviation is 0
        mean_d = np.mean(betti_control) - np.mean(betti_disease)
        if pooled_stdev == 0:
            observed_effect = np.nan  # Use NaN instead of 0
        else:
            observed_effect = mean_d / pooled_stdev
    
    elif effect_type == "ks":
        observed_effect = ks_2samp(betti_control, betti_disease).statistic

    else:
        raise ValueError("Invalid effect type. Choose 'cohen_d' or 'ks'.")

    # Combine data and permute
    combined = np.concatenate([betti_control, betti_disease])
    permuted_effects = []

    for _ in range(num_permutations):
        np.random.shuffle(combined)
        perm_control = combined[:len(betti_control)]
        perm_disease = combined[len(betti_control):]

        if effect_type == "cohen_d":
            pooled_stdev_perm = pooled_std(perm_control, perm_disease)
            d_perm = np.mean(perm_control) - np.mean(perm_disease)
            perm_effect = d_perm / pooled_stdev_perm
        else:
            perm_effect = ks_2samp(perm_control, perm_disease).statistic
            mean_d = np.nan

        permuted_effects.append(perm_effect)

    permuted_effects = np.array(permuted_effects)

    # Compute empirical p-value
    empirical_p = np.mean(np.abs(permuted_effects) >= np.abs(observed_effect))

    return observed_effect, mean_d, empirical_p

# -------------------
# Homology
# -------------------

def betti_to_pairs_single(betti_numbers, filtration_values):
    persistent_pairs = []  # To store (birth, death) pairs
    births = []  # To store birth times

    for i, current_betti in enumerate(betti_numbers):
        if i == 0:
            # At the first filtration value, initialize all births
            births = [filtration_values[i]] * current_betti
        else:
            previous_betti = betti_numbers[i - 1]
            if current_betti > previous_betti:
                # New features born
                births.extend([filtration_values[i]] * (current_betti - previous_betti))
            elif current_betti < previous_betti:
                # Some features died
                for _ in range(previous_betti - current_betti):
                    birth_time = births.pop(0)
                    persistent_pairs.append((birth_time, filtration_values[i]))
    
    # Handle features that never die (i.e., their death is at infinity)
    for birth_time in births:
        persistent_pairs.append((birth_time, float('inf')))
    
    return persistent_pairs

def find_change_and_stabilization(betti, filtration):
    first_change = None
    stabilization_value = None

    for i in range(1, len(betti)):
        if first_change is None and betti[i] != betti[i-1]:
            first_change = filtration[i]
        
        # If stabilization happens immediately after the first change
        if first_change is not None and all(betti[i:] == [betti[i]] * (len(betti) - i)):
            if first_change == filtration[i]:
                stabilization_value = filtration[i]
            else:
                stabilization_value = filtration[i-1]
            break

    return first_change, stabilization_value

def plot_betti_numbers(path_id, betti_control=None, betti_disease=None, dim=0, dist_type='', save_plot=False):    
    # Determine how many subplots are needed based on provided data
    n_plots = 0
    if betti_control is not None:
        n_plots += 1
    if betti_disease is not None:
        n_plots += 1

    if n_plots == 0:
        print("Both betti_control and betti_disease are None. Nothing to plot.")
        return

    # Set up the figure with the appropriate number of rows
    fig, axs = plt.subplots(n_plots, 1, figsize=(12, 4 * n_plots), sharex=True)

    # Ensure axs is always iterable, even if there's only one plot
    if n_plots == 1:
        axs = [axs]

    # Helper function to plot a single Betti number set
    def plot_single_betti(ax, betti_values, title):
        # Scale x-axis positions to match desired range
        x_range = len(betti_values) / 100  # Assuming 100 steps per filtration unit
        x_positions = np.linspace(0, x_range, len(betti_values))
        bar_width = x_range / len(betti_values)  # Adjust bar width to fit the scaled x-axis
        ax.bar(x_positions, betti_values, color='lightblue', width=bar_width, align='center')
        
        # Finding first change and stabilization points (assuming they use raw indices)
        first_change, stabilization = find_change_and_stabilization(betti_values, x_positions)
        
        # Annotate first change point if found
        if first_change is not None:
            ax.axvline(x=first_change, color='green', linestyle='--')
            ax.annotate(f'{first_change:.2f}', 
                        xy=(first_change, max(betti_values)), 
                        xytext=(first_change - 0.1, max(betti_values) * 0.9),
                        color='green')
        
        # Annotate stabilization point if found
        if stabilization is not None:
            ax.axvline(x=stabilization, color='red', linestyle='--')
            ax.annotate(f'{stabilization:.2f}', 
                        xy=(stabilization, max(betti_values)), 
                        xytext=(stabilization + 0.025, max(betti_values) * 0.9),
                        color='red')

        ax.set_title(title)
        ax.set_ylabel('Betti Count')
        return x_range

    # Plot Betti numbers for Control if provided
    x_range_control = x_range_disease = 0
    if betti_control is not None:
        x_range_control = plot_single_betti(axs[0 if betti_disease is None else 0], betti_control, f'Betti Numbers for {path_id} (Control - Dim {dim})')
    
    # Plot Betti numbers for Disease if provided
    if betti_disease is not None:
        x_range_disease = plot_single_betti(axs[0 if betti_control is None else 1], betti_disease, f'Betti Numbers for {path_id} (Disease - Dim {dim})')

    # Set shared x-axis ticks and labels
    overall_x_range = max(x_range_control, x_range_disease)  # Use the larger range for consistent scaling
    num_ticks = 10  # Adjust number of ticks for readability
    tick_positions = np.linspace(0, overall_x_range, num_ticks)
    axs[-1].set_xticks(tick_positions)
    axs[-1].set_xticklabels([f'{t:.2f}' for t in tick_positions])
    axs[-1].set_xlabel(f'Filtration Value (0–{overall_x_range:.2f})')

    plt.tight_layout()

    # Save the plot if requested
    if save_plot:
        plt.savefig(f'Plots/betti_numbers_plots_{path_id}_{dist_type}_dim{dim}.png')
        print(f'Plot for {path_id} ({dist_type}-{dim}) saved successfully')
        plt.close()
    else:
        plt.show()

def plot_persistent_diagram_single(ax, persistent_pairs, path_id, label_type, dim=None):
    """
    Function to plot a single persistence diagram on a given axis.
    
    Parameters:
    - ax: The axis object to plot on.
    - pers_pairs: List of persistent pairs [(birth, death), (birth, death), ...]
    - path_id: Identifier for the path.
    - label_type: Control or Disease label.
    """
    all_deaths = []  # Get all non-infinity deaths across all dimensions
    
    # Iterate over pairs and extract non-infinity death values
    for birth, death in persistent_pairs:
        if death != float('inf'):
            all_deaths.append(death)

    # Find the maximum non-infinite death value, or set a large default if all are infinite
    if all_deaths:
        max_death = max(all_deaths)
    else:
        max_death = 10  # Set an arbitrary large value as a substitute for infinity

    # Plot each dimension's birth and death
    for birth, death in persistent_pairs:
        death = death if death != float('inf') else max_death * 2

        # Scatter plot for each dimension
        ax.scatter(birth, death, s=50, color='red')

        # Plot line for infinity replacement (∞)
        if death == max_death * 2:
            inf_replacement = max_death * 2  # Replace infinity with large value
            max_axes = max_death * 2.1  # Extend axis limit for visual clarity
                
            # Plot horizontal dashed line for infinity
            ax.plot([-0.1 * max_axes, max_axes], [inf_replacement, inf_replacement], linestyle='--', color='grey', linewidth=1.5)
                
            # Annotate infinity (∞)
            mid_x = max_axes / 2
            annotation_y = inf_replacement * 0.96  # Adjust vertical position
            ax.annotate('∞', xy=(mid_x, annotation_y), xytext=(mid_x, annotation_y),
                        fontsize=14, color='k', ha='center', va='bottom')

    # Set limits to ensure the axes are equal
    max_value = max_death * 2.1
    ax.set_xlim((-0.025 * max_value, max_value))
    ax.set_ylim((-0.025 * max_value, max_value))

    # Add a diagonal line indicating birth = death
    ax.plot([0, max_value], [0, max_value], 'k--')

    # Set labels and title
    ax.set_xlabel('Birth')
    ax.set_ylabel('Death')
    if dim is not None:
        ax.set_title(f'{path_id} Persistence Diagram ({label_type.capitalize()} - Dim {dim})')
    else:
        ax.set_title(f'{path_id} Persistence Diagram ({label_type.capitalize()})')

def plot_persistent_diagrams(path_id, pers_pairs_control=None, pers_pairs_disease=None, dim=None, distance_type='absolute', save_plot=False):
    """
    Function to plot two persistence diagrams side by side for control and disease.
    
    Parameters:
    - path_id: The path ID being plotted.
    - pers_pairs_control: Persistent pairs for the control group.
    - pers_pairs_disease: Persistent pairs for the disease group.
    - save_plot: Boolean indicating whether to save the plot as a file.
    - dim: Dimension of the PH (int).
    """

    n_plots = 0
    if pers_pairs_control is not None:
        n_plots += 1
    if pers_pairs_disease is not None:
        n_plots += 1

    if n_plots == 0:
        print("Both pers_pairs_control and pers_pairs_disease are None. Nothing to plot.")
        return

    # Create figure with 1 row and n_plots columns
    fig, axs = plt.subplots(1, n_plots, figsize=(10, 2.5 * n_plots), sharex=True)

    # Ensure axs is always iterable, even if there's only one plot
    if n_plots == 1:
        axs = [axs]
        
    # Plot for Control
    if pers_pairs_control is not None:
        plot_persistent_diagram_single(axs[0], pers_pairs_control, path_id, 'control', dim)

    # Plot for Disease
    if pers_pairs_disease is not None:
        if n_plots == 1:
            plot_persistent_diagram_single(axs[0], pers_pairs_disease, path_id, 'disease', dim)
        else:
            plot_persistent_diagram_single(axs[1], pers_pairs_disease, path_id, 'disease', dim)

    # Adjust layout
    plt.tight_layout()

    # Save plot if requested
    if save_plot and dim is not None:
        plt.savefig(f'Plots/{distance_type.capitalize()}/Dim{dim}/persistent_diagrams_for_{path_id}_{dim}.png')
        print(f'Persistence diagrams for {path_id} (dim {dim}) saved successfully')
        plt.close()
    else:
        # Show plot
        plt.show()

def plot_persistent_barcode_single(ax, pers_pairs, path_id, label_type, dim=None, uniform_length=1):
    """
    Function to plot a single persistence barcode on a given axis.
    
    Parameters:
    - ax: The axis object to plot on.
    - pers_pairs: List of persistent pairs [(birth, death), (birth, death), ...]
    - path_id: Identifier for the path.
    - label_type: Control or Disease label.
    - dim: Dimension of the PH (optional).
    """
    # Number of persistent pairs
    max_pairs = len(pers_pairs)
    multiplier = max_pairs / 15

    # Get all non-infinity deaths to calculate max_death
    all_deaths = [death for _, death in pers_pairs if death != float('inf')]

    if all_deaths:
        max_death = uniform_length #max(all_deaths, default=0)
    else:
        max_death = 1

    # Plot the pairs
    for idx, (birth, death) in enumerate(sorted(pers_pairs)):
        if death == float('inf'):
            # Replace infinity with a large value and add an arrow
            inf_replacement = max_death * 2
            death = inf_replacement
            ax.plot([birth, death], [idx, idx], color=f"C0", lw=2)
            ax.annotate('', xy=(death, idx), xytext=(death - 0.1, idx),
                        arrowprops=dict(arrowstyle="->", color=f"C0"))
        else:
            ax.plot([birth, death], [idx, idx], color=f"C0", lw=2)
    
    # Adjust plot settings
    ax.set_ylim([-0.5, max_pairs - 0.5])
    ax.set_xlabel("Filtration Value")
    
    # Set y-axis ticks spaced by the given step size (y_step)
    y_ticks = np.arange(0, max_pairs, 10)
    ax.set_yticks(y_ticks)
    
    if dim is not None:
        ax.set_title(f'{path_id} Barcode ({label_type.capitalize()} - Dimension {dim})')
    else:
        ax.set_title(f'{path_id} Barcode ({label_type.capitalize()})')

def plot_persistent_barcode(path_id, pers_pairs_control=None, pers_pairs_disease=None, dim=None, distance_type='', save_plot=False):
    """
    Function to plot two persistence barcodes side by side for control and disease.
    
    Parameters:
    - path_id: The path ID being plotted.
    - pers_pairs_control: Persistent pairs for the control group.
    - pers_pairs_disease: Persistent pairs for the disease group.
    - save_plot: Boolean indicating whether to save the plot as a file.
    - dim: Dimension of the PH (int).
    """

    n_plots = 0
    if pers_pairs_control is not None:
        n_plots += 1
    if pers_pairs_disease is not None:
        n_plots += 1

    if n_plots == 0:
        print("Both pers_pairs_control and pers_pairs_disease are None. Nothing to plot.")
        return

    # Create figure with 1 row and n_plots columns
    fig, axs = plt.subplots(n_plots, 1, figsize=(10, 5 * n_plots), sharex=True)

    # Ensure axs is always iterable, even if there's only one plot
    if n_plots == 1:
        axs = [axs]

    # Determine uniform length for x-axis
    uniform_length = 1  # Default value in case of no finite deaths
    if pers_pairs_control is not None:
        all_deaths_control = [death for _, death in pers_pairs_control if death != float('inf')]
        if all_deaths_control:
            uniform_length = max(uniform_length, max(all_deaths_control))
    if pers_pairs_disease is not None:
        all_deaths_disease = [death for _, death in pers_pairs_disease if death != float('inf')]
        if all_deaths_disease:
            uniform_length = max(uniform_length, max(all_deaths_disease))
        
    # Plot for Control
    if pers_pairs_control is not None:
        plot_persistent_barcode_single(axs[0], pers_pairs_control, path_id, 'control', dim, uniform_length)

    # Plot for Disease
    if pers_pairs_disease is not None:
        if n_plots == 1:
            plot_persistent_barcode_single(axs[0], pers_pairs_disease, path_id, 'disease', dim, uniform_length)
        else:
            plot_persistent_barcode_single(axs[1], pers_pairs_disease, path_id, 'disease', dim, uniform_length)

    # Adjust layout
    plt.tight_layout()

    # Save plot if requested
    if save_plot and dim is not None:
        plt.savefig(f'Plots/{distance_type.capitalize()}/Dim{dim}/persistent_barcode_for_{path_id}_dim_{dim}.png')
        #plt.savefig(f'Plots/Euclidean/Dim0/persistent_barcode_for_{path_id}_dim_{dim}.png')
        print(f'Persistent barcode for {path_id} (dim {dim}) saved successfully')
        plt.close()
    else:
        # Show plot
        plt.show()

