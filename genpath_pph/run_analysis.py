"""
run_analysis.py
================

This script performs statistical analyses on persistence-based pathway homology (PPH) results. 
It provides methods to quantify and test topological differences between control and disease 
groups based on persistence landscapes and Betti number distributions.

Main functionalities:
---------------------

1. pooled_std:
   - Computes the pooled standard deviation between two groups, used in effect size calculations.

2. permutation_test_effect_size:
   - Estimates the effect size (Cohen’s d or Kolmogorov–Smirnov statistic) between 
     Betti number distributions of control and disease groups.
   - Performs permutation testing to obtain empirical p-values.

3. permutation_test:
   - Conducts a permutation test on persistence landscapes using various norms (L1, L2, sup)
     to assess global topological differences.

4. global_difference_analysis:
   - Aggregates results of the permutation tests across multiple norms.
   - Applies Benjamini–Hochberg correction for multiple testing and reports 
     significant global topological differences.

5. pathway_level_difference_analysis:
   - Performs pathway-level testing using effect size and KS statistics on Betti numbers.
   - Applies permutation-based empirical p-value estimation and multiple testing correction.
   - Saves the resulting statistics to CSV files in the 'Results/' directory.

Dependencies:
-------------
time, numpy, pandas, random, pickle, scipy.stats, statsmodels.stats.multitest

Output:
-------
- Global difference summary (printed to console)
- Pathway-level test results saved as:
  'Results/Effect_size_KS_test_results_dimension_{top_dim}_{dist_type}_new.csv'
"""

import time                     # Used to track or measure execution time of operations
import numpy as np              # Numerical computations, array and matrix handling
import pandas as pd             # For data handling and manipulation
import random                   # Random sampling for permutation testing
from scipy.stats import ks_2samp            # For performing the KS test
from statsmodels.stats.multitest import multipletests  # Multiple testing correction (e.g., Benjamini–Hochberg)
from persim.landscapes import average_approx, snap_pl  # Averaging and aligning (snapping) landscapes

def permutation_test(cpls, dpls, num_perms=1000, norm=1):
    np.random.seed(0)

    num_runs = len(cpls)
    combined = cpls + dpls

    [c_sn, d_sn] = snap_pl([average_approx(cpls), average_approx(dpls)])
    true_diff = (d_sn - c_sn)
    true_diff_val = (
        true_diff.p_norm(p=norm) if norm in [1, 2] else true_diff.sup_norm()
    )

    sig_count = 0
    for i in range(num_perms):
        # Print progress every 10 permutations
        if (i + 1) % 100 == 0:
            print(f"Permutation {i + 1}/{num_perms} completed (for {norm}-norm)")
        
        A_indices = random.sample(range(2 * num_runs), num_runs)
        B_indices = [_ for _ in range(2*num_runs) if _ not in A_indices]
        A = [combined[i] for i in A_indices]
        B = [combined[j] for j in B_indices]
        [A_sn, B_sn] = snap_pl([average_approx(A), average_approx(B)])
        
        random_diff = (B_sn - A_sn)
        diff_val = (
            random_diff.p_norm(p=norm) if norm in [1, 2] else random_diff.sup_norm()
        )
        if diff_val >= true_diff_val:
            sig_count += 1
    return true_diff_val, (sig_count / num_perms)

def pooled_std(control, disease):
    n1, n2 = len(control), len(disease)
    var1, var2 = np.var(control, ddof=1), np.var(disease, ddof=1)
    return np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    
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
    
def global_difference_analysis(cpls, dpls, num_perms=1000):
    # Snap & compute averages
    cpls_snap, dpls_snap = [], []
    for cpl, dpl in zip(cpls, dpls):
        try:
            [c_sn, d_sn] = snap_pl([cpl, dpl])
            cpls_snap.append(c_sn)
            dpls_snap.append(d_sn)
        except Exception:
            pass
            
    # Run permutation test for 3 norms
    norms = [1, 2, 'sup']
    results_raw = [permutation_test(cpls_snap, dpls_snap, num_perms=num_perms, norm=n) for n in norms]
    true_diff_vals, pvals = zip(*results_raw)
    
    # BH correction
    reject, pvals_adj, _, _ = multipletests(pvals, method='fdr_bh')
    results = dict(zip(norms, zip(pvals, pvals_adj, reject)))

    # Display the global difference test results
    print("\nGlobal difference analysis:")
    for i, norm in enumerate(norms):
        p_raw, p_adj, sig = results[norm]
        diff_val = true_diff_vals[i]
        print(f"Norm {norm}: measured difference={diff_val:.6f}, p={p_raw:.6f}, adj p={p_adj:.6f}, Significant={sig}")
        if p_adj < 0.05:
            print(f"{norm}-norm shows significant difference between control and disease groups\n")
        else:
            print(f"{norm}-norm did not show significant difference between control and disease groups\n")
    
    return results

def pathway_level_difference_analysis(path_ids, dim0_betti_numbers_dict, dim1_betti_numbers_dict, top_dim=0, num_perm=5000, seed = 0):
    start_time = time.time()
    
    dist_type = '1-abs-correlation' # Distance function used in GenPath-PPH

    # Check top_dim and dist_type to select the appropriate dictionary
    if top_dim == 0:
        dim_betti_numbers_dict = dim0_betti_numbers_dict
    elif top_dim == 1:
        dim_betti_numbers_dict = dim1_betti_numbers_dict
    else:
        raise ValueError(f"Invalid entry for top_dim ({top_dim}).")

    # Initialize lists for columns
    dim_results_data = {
        'path_id': [],
        'mean_diff_observed': [],
        'es_observed': [],
        'es_raw_pvalue': [],
        'es_corrected_pvalue': [],
        'ks_observed': [],
        'ks_raw_pvalue': [],
        'ks_corrected_pvalue': []
        }

    # Loop over each path_id to perform KS test for full, cut_both_ends, and cut_tail_end
    counter = 0
    for path_id in path_ids:
        counter += 1
        print(f"Processing pathway {counter}/{len(path_ids)}: {path_id}...")

        dim_results_data['path_id'].append(path_id)
    
        # Full Betti numbers
        betti_control_full = dim_betti_numbers_dict[path_id]['control']
        betti_disease_full = dim_betti_numbers_dict[path_id]['disease']
        observed_diff_es, mean_d, raw_p_es = permutation_test_effect_size(betti_control_full, betti_disease_full, num_permutations=num_perm, effect_type="cohen_d", seed=seed)
        observed_diff_ks, _, raw_p_ks = permutation_test_effect_size(betti_control_full, betti_disease_full, num_permutations=num_perm, effect_type="ks", seed=seed)
        dim_results_data['mean_diff_observed'].append(mean_d)
        dim_results_data['es_observed'].append(observed_diff_es)
        dim_results_data['es_raw_pvalue'].append(raw_p_es)
        dim_results_data['ks_observed'].append(observed_diff_ks)
        dim_results_data['ks_raw_pvalue'].append(raw_p_ks)

    # Apply Benjamini-Hochberg correction on raw p-values for each group
    _, es_corrected_pvalues, _, _ = multipletests(dim_results_data['es_raw_pvalue'], method='fdr_bh')
    _, ks_corrected_pvalues, _, _ = multipletests(dim_results_data['ks_raw_pvalue'], method='fdr_bh')

    # Store the corrected p-values
    dim_results_data['es_corrected_pvalue'] = es_corrected_pvalues
    dim_results_data['ks_corrected_pvalue'] = ks_corrected_pvalues

    # Convert the results to a dataframe
    dim_df_results = pd.DataFrame(dim_results_data)

    # save the resulting dataframe
    dim_df_results.to_csv(f"Results/Effect_size_KS_test_results_dimension_{top_dim}_{dist_type}_new.csv", index=False)

    print(f"KS test run for {time.time() - start_time} secs.")
    
