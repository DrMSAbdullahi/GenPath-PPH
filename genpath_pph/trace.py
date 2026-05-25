"""
genpath_pph/trace.py
====================
Inverse mapping functions for GenPath-PPH persistent path homology results.

Given the Betti number series and edges-at-each-filtration output from
persistent_path_homology_from_digraph (core.py), these functions trace
which genes and gene-gene interactions are responsible for each
persistent topological feature.

Mathematical basis
------------------
β₀ — Connected components
    Detected via Union-Find on undirected projection of active edges,
    consistent with split_independent_component in core.py.

β₁ — Path homology cycles
    A path homology 1-cycle is a vector c in ker(∂₁) ∩ Ω₁ that is not
    in im(∂₂) ∩ Ω₁, where:
        Ω₁  = space of allowed 1-paths (active pathway edges)
        ∂₁  = boundary operator (1-paths → 0-paths)
        ∂₂  = boundary operator (2-paths → 1-paths)

    IMPORTANT: Path homology cycles are NOT simply directed graph loops.
    A diamond pattern g1→g2, g1→g3, g2→g4, g3→g4 IS a valid path
    homology 1-cycle (the formal sum g1→g2→g4 - g1→g3→g4 has zero
    boundary) even though it contains no directed loop.

    Implementation uses the extended space formulation from core.py:
        dim(im(∂₂) ∩ Ω₁) = dim(Ω₁) + rank(∂₂) - dim_An_Bn
    where dim_An_Bn = rank(vstack([eye * path_idx[2], ∂₂])).
    When this quantity is 0, all ker(∂₁) vectors are genuine cycles.
    When > 0, only the complement of im(∂₂) in ker(∂₁) gives true cycles.

Validation
----------
Tested on toy pathway with:
    - Directed loop g1→g2→g3→g4→g1 (detected at filtration=0.95)
    - Diamond g6→g7, g6→g8, g9→g7, g9→g8 (detected at filtration=0.39)
Both correctly identified as independent path homology 1-cycles.

Reference
---------
Abdullahi et al. (2025). GenPath-PPH. CSBJ, 27, 5348-5362.
https://doi.org/10.1016/j.csbj.2025.11.018

Chowdhury & Memoli (2018). Persistent path homology of directed networks.
SODA 2018.
"""

from __future__ import annotations

import warnings
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ══════════════════════════════════════════════════════════════════════════════
# Internal: edge validation
# ══════════════════════════════════════════════════════════════════════════════

def validate_edges_at_filtration(
    edges_at_each_filtration: list,
    adj_matrix: np.ndarray,
    filtration: np.ndarray,
    n_nodes: int,
) -> dict:
    """
    Sanity check on edges_at_each_filtration from core.py.

    Verifies four properties:
    1. Edge count is monotonically non-decreasing (edges only added)
    2. All node indices in valid range [0, n_nodes)
    3. All edges exist in the KEGG pathway adjacency matrix
    4. No self-loops (i == j)

    Parameters
    ----------
    edges_at_each_filtration : list
        Second return value of persistent_path_homology_from_digraph.
    adj_matrix : np.ndarray, shape (n_nodes, n_nodes)
        Binary directed adjacency matrix of the KEGG pathway graph.
    filtration : np.ndarray
        Filtration threshold values.
    n_nodes : int
        Number of genes in the pathway.

    Returns
    -------
    dict with keys:
        'monotone', 'valid_indices', 'pathway_compliant',
        'no_self_loops', 'edge_counts', 'violations'

    Example
    -------
        val = validate_edges_at_filtration(
            edges_filt, adj_matrix.values, filtration, n_nodes=len(gene_names)
        )
        # ✓ Edge validation passed — all checks OK.
    """
    violations = []
    edge_counts = [len(e) for e in edges_at_each_filtration]

    # Check 1 — monotone non-decreasing
    monotone = all(
        edge_counts[i] <= edge_counts[i+1]
        for i in range(len(edge_counts)-1)
    )
    if not monotone:
        drops = [i for i in range(len(edge_counts)-1)
                 if edge_counts[i] > edge_counts[i+1]]
        violations.append(
            f"Edge count decreases at steps: {drops}. "
            f"Edges must be added monotonically as filtration increases."
        )

    # Checks 2, 3, 4
    valid_indices     = True
    pathway_compliant = True
    no_self_loops     = True

    for t, edges in enumerate(edges_at_each_filtration):
        if len(edges) == 0:
            continue
        for row in np.array(edges, dtype=int):
            i, j = int(row[0]), int(row[1])
            if i < 0 or i >= n_nodes or j < 0 or j >= n_nodes:
                valid_indices = False
                violations.append(
                    f"Step {t}: index ({i},{j}) out of range [0,{n_nodes})."
                )
            if i == j:
                no_self_loops = False
                violations.append(f"Step {t}: self-loop at node {i}.")
            if (0 <= i < n_nodes and 0 <= j < n_nodes
                    and adj_matrix[i, j] == 0):
                pathway_compliant = False
                violations.append(
                    f"Step {t}: edge ({i}→{j}) not in pathway graph."
                )

    result = {
        "monotone":          monotone,
        "valid_indices":     valid_indices,
        "pathway_compliant": pathway_compliant,
        "no_self_loops":     no_self_loops,
        "edge_counts":       edge_counts,
        "violations":        violations,
    }

    if not violations:
        print("✓ Edge validation passed — all checks OK.")
    else:
        print(f"✗ {len(violations)} violation(s) found:")
        for v in violations[:10]:
            print(f"  - {v}")

    return result


# ══════════════════════════════════════════════════════════════════════════════
# Internal: β₀ — connected components via Union-Find
# ══════════════════════════════════════════════════════════════════════════════

def _union_find_components(
    edges: np.ndarray,
    n_nodes: int,
) -> Dict[int, List[int]]:
    """
    Find weakly connected components using Union-Find.

    Treats directed edges as undirected, consistent with
    split_independent_component in core.py (undirected DFS).
    """
    parent = list(range(n_nodes))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    if len(edges) > 0:
        for row in edges:
            union(int(row[0]), int(row[1]))

    components = defaultdict(list)
    for node in range(n_nodes):
        components[find(node)].append(node)

    return dict(components)


# ══════════════════════════════════════════════════════════════════════════════
# Internal: β₁ — path homology cycles via boundary operator null space
# ══════════════════════════════════════════════════════════════════════════════

def _get_ph_cycles_at_step(
    edges_arr: np.ndarray,
    n_nodes: int,
    gene_names: List[str],
    pph_instance,
    tol: float = 1e-8,
) -> List[Dict]:
    """
    Find all independent path homology 1-cycles at one filtration step.

    Uses utils_generate_allowed_paths and utils_unlimited_boundary_operator
    from core.py exactly as they are used in path_homology_for_connected_digraph.

    Mathematical approach
    ---------------------
    1. Compute ker(∂₁): null space of the boundary operator on 1-paths.
       Each vector in ker(∂₁) is a candidate cycle — a linear combination
       of allowed edges whose boundary (sum of endpoints) is zero.

    2. Compute dim(im(∂₂) ∩ Ω₁) using the extended space formula from core.py:
           dim_1 = len(Ω₁) + rank(∂₂) - rank(vstack([eye*path_idx[2], ∂₂]))
       This measures how many ker(∂₁) vectors are actually boundaries of
       2-paths and therefore NOT genuine homology classes.

    3. If dim_1 = 0: ALL ker(∂₁) vectors are genuine cycles (no projection).
       If dim_1 > 0: project out im(∂₂) from ker(∂₁) using QR decomposition
       on the correctly reindexed ∂₂ matrix.

    Parameters
    ----------
    edges_arr : np.ndarray, shape (n_edges, 2)
        Active edges at this filtration step (integer node indices).
    n_nodes : int
    gene_names : list of str
    pph_instance : GenPathHomology instance from core.py
    tol : float, default 1e-8

    Returns
    -------
    list of dicts, each with:
        'gene_members' : list of gene symbols in this cycle
        'edges'        : list of (gene_i, gene_j) edge tuples
        'paths'        : list of path strings (e.g. 'g1->g2')
        'coefficients' : dict mapping path string to coefficient
        'cycle_vector' : raw null space vector (one entry per 1-path)
        'n_genes'      : number of genes in this cycle
    """
    if len(edges_arr) == 0:
        return []

    # ── Step 1: generate allowed paths up to dimension 2 ──────────────────
    try:
        allowed_paths = pph_instance.utils_generate_allowed_paths(
            edges_arr, max_path=2
        )
    except Exception as e:
        warnings.warn(f"utils_generate_allowed_paths failed: {e}")
        return []

    one_paths = allowed_paths.get(1, [])
    if not one_paths:
        return []

    n_1paths = len(one_paths)

    # ── Step 2: compute boundary matrices ─────────────────────────────────
    try:
        boundary_matrices, ranks, path_idx = \
            pph_instance.utils_unlimited_boundary_operator(
                allowed_paths, max_path=2
            )
    except Exception as e:
        warnings.warn(f"utils_unlimited_boundary_operator failed: {e}")
        return []

    D1 = boundary_matrices.get(1)   # shape (n_1paths, n_considered_0paths)
    D2 = boundary_matrices.get(2)   # shape (n_2paths, n_considered_1paths)

    if D1 is None or D1.shape[0] == 0:
        return []

    # ── Step 3: find ker(∂₁) ──────────────────────────────────────────────
    # Vectors c (length n_1paths) such that D1 @ c = 0
    # = right null space of D1 = null space of D1.T
    D1_float = D1.astype(float)

    try:
        _, S1, Vt1 = np.linalg.svd(D1_float.T, full_matrices=True)
    except np.linalg.LinAlgError as e:
        warnings.warn(f"SVD of D1 failed: {e}")
        return []

    ker_basis = []
    for i in range(len(S1)):
        if S1[i] < tol:
            ker_basis.append(Vt1[i])
    for i in range(len(S1), Vt1.shape[0]):
        ker_basis.append(Vt1[i])

    if not ker_basis:
        return []

    # ── Step 4: compute dim(im(∂₂) ∩ Ω₁) using core.py formula ──────────
    # From path_homology_for_connected_digraph in core.py:
    #   dim_An_Bn = rank(vstack([eye * path_idx[2], D2]))
    #   dim_1 = len(Ω₁) + rank(D2) - dim_An_Bn
    # dim_1 measures how many ker vectors are boundaries (not genuine cycles)

    dim_1 = 0
    im_basis_9dim = None   # orthonormal basis for im(∂₂) in Ω₁ if needed

    if D2 is not None and D2.shape[0] > 0:
        flags_2  = path_idx.get(2, [1] * D2.shape[1])
        D2_float = D2.astype(float)

        # Replicate core.py formula exactly
        eye_part = np.eye(len(flags_2)) * flags_2
        stacked  = np.vstack([eye_part, D2_float])
        dim_An_Bn = int(np.linalg.matrix_rank(stacked))
        rank_D2   = ranks.get(2, int(np.linalg.matrix_rank(D2_float)))

        dim_1 = n_1paths + rank_D2 - dim_An_Bn

    # ── Step 5: get true cycle generators ─────────────────────────────────
    if dim_1 == 0:
        # No elements of ker(∂₁) are in im(∂₂) — all ker vectors genuine
        true_cycles = ker_basis

    else:
        # Project out im(∂₂) from ker(∂₁)
        # Build ∂₂ in the correct 9-dim Ω₁ basis using reindex map
        # (∂₂ columns use sorted-path ordering, Ω₁ uses allowed_paths ordering)

        # Reconstruct the considered column ordering of ∂₂
        two_paths  = allowed_paths.get(2, [])
        flags_2    = path_idx.get(2, [1] * D2.shape[1])

        boundary_path_names = []
        for tpath in two_paths:
            nodes = tpath.split("->")
            for i_kill in range(len(nodes)):
                remaining = [n for k, n in enumerate(nodes) if k != i_kill]
                boundary_path_names.append("->".join(remaining))

        considered_cols = sorted(set(boundary_path_names + one_paths))

        # Reindex: D2 column index → Ω₁ index (position in one_paths)
        reindex = {}
        for j, path in enumerate(considered_cols):
            if j < len(flags_2) and flags_2[j] == 1 and path in one_paths:
                reindex[j] = one_paths.index(path)

        # Build ∂₂ restricted to Ω₁ columns in Ω₁ ordering
        D2_omega1 = np.zeros((D2.shape[0], n_1paths))
        for d2_col, omega1_idx in reindex.items():
            D2_omega1[:, omega1_idx] = D2_float[:, d2_col]

        # Orthonormal basis for im(∂₂) in Ω₁ via SVD
        try:
            Uim, Sim, Vtim = np.linalg.svd(D2_omega1, full_matrices=False)
            im_basis = Vtim[Sim > tol]   # rows = orthonormal basis vectors
        except np.linalg.LinAlgError:
            im_basis = np.empty((0, n_1paths))

        # Project ker vectors onto complement of im(∂₂)
        true_cycles = []
        for vec in ker_basis:
            if len(im_basis) > 0:
                proj      = sum(np.dot(vec, b) * b for b in im_basis)
                remainder = vec - proj
            else:
                remainder = vec
            if np.linalg.norm(remainder) > tol:
                true_cycles.append(remainder)

    # ── Step 6: map cycle vectors to gene-level information ───────────────
    cycles_result = []
    seen = set()

    for vec in true_cycles:
        if len(vec) != n_1paths:
            continue

        active_idx = [i for i, v in enumerate(vec) if abs(v) > tol]
        if not active_idx:
            continue

        # Deduplicate by active edge set
        key = frozenset(active_idx)
        if key in seen:
            continue
        seen.add(key)

        cycle_edges  = []
        cycle_genes  = set()
        cycle_paths  = []
        coefficients = {}

        for idx in active_idx:
            path_str  = one_paths[idx]        # e.g. "3->7"
            node_strs = path_str.split("->")

            try:
                node_ids = [int(n) for n in node_strs]
            except ValueError:
                continue

            if not all(n < len(gene_names) for n in node_ids):
                continue

            gene_path  = [gene_names[n] for n in node_ids]
            path_label = "->".join(gene_path)
            coeff      = round(float(vec[idx]), 4)

            cycle_paths.append(path_label)
            coefficients[path_label] = coeff
            cycle_genes.update(gene_path)

            if len(gene_path) == 2:
                cycle_edges.append((gene_path[0], gene_path[1]))

        if not cycle_genes:
            continue

        cycles_result.append({
            "gene_members": sorted(cycle_genes),
            "edges":        cycle_edges,
            "paths":        cycle_paths,
            "coefficients": coefficients,
            "cycle_vector": vec.tolist(),
            "n_genes":      len(cycle_genes),
        })

    return cycles_result


# ══════════════════════════════════════════════════════════════════════════════
# Main function: trace_persistent_features
# ══════════════════════════════════════════════════════════════════════════════

def trace_persistent_features(
    betti_series: list,
    edges_at_each_filtration: list,
    n_nodes: int,
    gene_names: List[str],
    pph_instance=None,
    target_dimension: int = 0,
    filtration: Optional[np.ndarray] = None,
    min_component_size: int = 2,
    verbose: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Trace which genes are involved in each persistent topological feature.

    For β₀: identifies gene members of each weakly connected component
    at every filtration step where the Betti number changes. A change
    means components form (edges first connect nodes) or merge.

    For β₁: identifies genes and edges forming each path homology cycle
    using the null space of the boundary operator ∂₁, with im(∂₂)
    projected out using the extended space formulation from core.py.
    Both directed loops AND diamond patterns are correctly identified.

    Parameters
    ----------
    betti_series : list of int
        Betti series from persistent_path_homology_from_digraph.
    edges_at_each_filtration : list of np.ndarray
        Edge sets at each filtration step (second return value of
        persistent_path_homology_from_digraph). Each entry shape
        (n_active_edges, 2) with integer node indices.
    n_nodes : int
        Number of genes in the pathway.
    gene_names : list of str
        gene_names[i] = gene symbol for node index i.
        Must match row order of expression matrix passed to PPH.
    pph_instance : GenPathHomology, required for target_dimension=1
        An instance of GenPathHomology from core.py.
    target_dimension : int, default 0
        0 for component tracing (β₀), 1 for cycle tracing (β₁).
    filtration : np.ndarray, optional
        Filtration threshold values for labelling output.
    min_component_size : int, default 2
        Minimum genes in a component to report (β₀ only).
        Isolated single-gene nodes excluded by default.
    verbose : bool, default False
        Print progress at each filtration step where β changes.

    Returns
    -------
    features_df : pd.DataFrame
        One row per feature per filtration step where Betti changes.
        Columns:
            step, filtration_value, betti, change,
            n_features, feature_id, feature_type,
            gene_members, n_genes,
            edges_in_feature, n_edges,
            paths (β₁ only), coefficients (β₁ only)

    summary_df : pd.DataFrame
        One row per filtration step where Betti changes.
        Columns:
            step, filtration_value, betti, change,
            n_features, all_genes_involved, n_unique_genes

    Raises
    ------
    ValueError
        If gene_names length != n_nodes, or series lengths mismatch.
    RuntimeError
        If target_dimension=1 and pph_instance is None.

    Example
    -------
        from genpath_pph import GenPathHomology, load_toy_data
        from genpath_pph.trace import trace_persistent_features

        df, edges, filtration = load_toy_data()
        gene_names = df.index.tolist()
        X = df.values
        pph = GenPathHomology()

        # β₀ — connected components
        betti_0, edges_filt = pph.persistent_path_homology_from_digraph(
            X, edges, target_dimension=0, filtration=filtration
        )
        feat_df, summ_df = trace_persistent_features(
            betti_0, edges_filt, len(gene_names), gene_names,
            target_dimension=0, filtration=filtration
        )

        # β₁ — path homology cycles (pass pph_instance)
        betti_1, edges_filt_1 = pph.persistent_path_homology_from_digraph(
            X, edges, target_dimension=1, filtration=filtration
        )
        feat_df_1, summ_df_1 = trace_persistent_features(
            betti_1, edges_filt_1, len(gene_names), gene_names,
            pph_instance=pph,
            target_dimension=1, filtration=filtration
        )
        print(summ_df_1)
    """
    if len(gene_names) != n_nodes:
        raise ValueError(
            f"gene_names length ({len(gene_names)}) must equal "
            f"n_nodes ({n_nodes})."
        )
    if len(betti_series) != len(edges_at_each_filtration):
        raise ValueError(
            f"betti_series and edges_at_each_filtration must have "
            f"the same length."
        )
    if target_dimension == 1 and pph_instance is None:
        raise RuntimeError(
            "pph_instance (GenPathHomology) is required for β₁ tracing. "
            "Pass pph_instance=GenPathHomology()."
        )

    if filtration is None:
        filtration = np.linspace(0, 1, len(betti_series))

    feature_rows = []
    summary_rows = []
    prev_betti   = int(betti_series[0])

    for t in range(len(betti_series)):
        curr_betti = int(betti_series[t])

        # Process only steps where Betti changes (always process step 0)
        if curr_betti == prev_betti and t > 0:
            prev_betti = curr_betti
            continue

        fval   = round(float(filtration[t]), 6)
        change = curr_betti - prev_betti if t > 0 else curr_betti

        edges_raw = edges_at_each_filtration[t]
        edges_arr = (np.array(edges_raw, dtype=int)
                     if len(edges_raw) > 0
                     else np.empty((0, 2), dtype=int))

        if verbose:
            print(f"  Step {t:3d} (filtration={fval:.4f}): "
                  f"β={curr_betti}, change={change:+d}, "
                  f"active edges={len(edges_arr)}")

        all_genes_this_step = []
        n_features_this_step = 0

        # ── β₀: connected components ──────────────────────────────────────
        if target_dimension == 0:
            components = _union_find_components(edges_arr, n_nodes)
            comp_list  = sorted(components.values(), key=len, reverse=True)
            comp_list  = [c for c in comp_list if len(c) >= min_component_size]

            n_features_this_step = len(comp_list)

            for feat_id, component in enumerate(comp_list):
                gene_members = [gene_names[n] for n in component]
                comp_set     = set(component)

                # Edges within this component
                within_edges = (edges_arr[np.array([
                    int(r[0]) in comp_set and int(r[1]) in comp_set
                    for r in edges_arr
                ])] if len(edges_arr) > 0 else [])

                edge_pairs = [
                    (gene_names[int(r[0])], gene_names[int(r[1])])
                    for r in within_edges
                ]
                all_genes_this_step.extend(gene_members)

                feature_rows.append({
                    "step":             t,
                    "filtration_value": fval,
                    "betti":            curr_betti,
                    "change":           change,
                    "n_features":       n_features_this_step,
                    "feature_id":       feat_id,
                    "feature_type":     "component",
                    "gene_members":     gene_members,
                    "n_genes":          len(gene_members),
                    "edges_in_feature": edge_pairs,
                    "n_edges":          len(edge_pairs),
                    "paths":            None,
                    "coefficients":     None,
                })

        # ── β₁: path homology cycles ──────────────────────────────────────
        elif target_dimension == 1:
            cycles = _get_ph_cycles_at_step(
                edges_arr, n_nodes, gene_names, pph_instance
            )
            n_features_this_step = len(cycles)

            if not cycles:
                feature_rows.append({
                    "step":             t,
                    "filtration_value": fval,
                    "betti":            curr_betti,
                    "change":           change,
                    "n_features":       0,
                    "feature_id":       None,
                    "feature_type":     "ph_cycle",
                    "gene_members":     [],
                    "n_genes":          0,
                    "edges_in_feature": [],
                    "n_edges":          0,
                    "paths":            [],
                    "coefficients":     {},
                })
            else:
                for feat_id, cycle in enumerate(cycles):
                    all_genes_this_step.extend(cycle["gene_members"])
                    feature_rows.append({
                        "step":             t,
                        "filtration_value": fval,
                        "betti":            curr_betti,
                        "change":           change,
                        "n_features":       n_features_this_step,
                        "feature_id":       feat_id,
                        "feature_type":     "ph_cycle",
                        "gene_members":     cycle["gene_members"],
                        "n_genes":          cycle["n_genes"],
                        "edges_in_feature": cycle["edges"],
                        "n_edges":          len(cycle["edges"]),
                        "paths":            cycle["paths"],
                        "coefficients":     cycle["coefficients"],
                    })

        # ── Summary Row ─────────────────────────────────────────────────
        summary_rows.append({
            "step":               t,
            "filtration_value":   fval,
            "betti":              curr_betti,
            "change":             change,
            "n_features":         n_features_this_step,
            "all_genes_involved": sorted(set(all_genes_this_step)),
            "n_unique_genes":     len(set(all_genes_this_step)),
        })

        prev_betti = curr_betti

    features_df = (pd.DataFrame(feature_rows)
                   if feature_rows else pd.DataFrame())
    summary_df  = (pd.DataFrame(summary_rows)
                   if summary_rows  else pd.DataFrame())

    return features_df, summary_df


# ══════════════════════════════════════════════════════════════════════════════
# Gene frequency summary
# ══════════════════════════════════════════════════════════════════════════════

def summarise_gene_frequency(
    features_df: pd.DataFrame,
    top_n: int = 20,
) -> pd.DataFrame:
    """
    Rank genes by how frequently they appear across persistent features.

    Genes appearing in many filtration steps are topologically persistent —
    they remain involved in components or cycles across a wide range of
    correlation thresholds, indicating structural importance in the pathway.

    Parameters
    ----------
    features_df : pd.DataFrame
        First return value of trace_persistent_features.
    top_n : int, default 20

    Returns
    -------
    pd.DataFrame columns:
        gene, n_steps, n_features, first_step, last_step, persistence

    Example
    -------
        freq = summarise_gene_frequency(feat_df, top_n=10)
        print(freq)
    """
    if features_df.empty:
        return pd.DataFrame()

    gene_stats = defaultdict(lambda: {
        "steps": set(), "n_features": 0, "all_steps": []
    })

    for _, row in features_df.iterrows():
        members = row.get("gene_members") or []
        for gene in members:
            gene_stats[gene]["steps"].add(row["step"])
            gene_stats[gene]["n_features"] += 1
            gene_stats[gene]["all_steps"].append(row["step"])

    records = []
    for gene, s in gene_stats.items():
        steps = s["all_steps"]
        records.append({
            "gene":        gene,
            "n_steps":     len(s["steps"]),
            "n_features":  s["n_features"],
            "first_step":  min(steps),
            "last_step":   max(steps),
            "persistence": max(steps) - min(steps),
        })

    return (pd.DataFrame(records)
            .sort_values("n_steps", ascending=False)
            .head(top_n)
            .reset_index(drop=True))


# ══════════════════════════════════════════════════════════════════════════════
# Compare disease vs control features
# ══════════════════════════════════════════════════════════════════════════════

def compare_group_features(
    betti_disease: list,
    edges_disease: list,
    betti_control: list,
    edges_control: list,
    n_nodes: int,
    gene_names: List[str],
    pph_instance=None,
    target_dimension: int = 0,
    filtration: Optional[np.ndarray] = None,
    top_n: int = 20,
) -> Dict:
    """
    Compare genes involved in persistent features between disease and control.

    Runs trace_persistent_features on both groups and identifies:
    - Disease-specific genes (in disease features, not in control)
    - Control-specific genes (in control features, not in disease)
    - Shared genes (in features of both groups)

    Parameters
    ----------
    betti_disease, betti_control : list — Betti series for each group
    edges_disease, edges_control : list — edges_at_each_filtration for each
    n_nodes : int
    gene_names : list of str
    pph_instance : GenPathHomology, required for target_dimension=1
    target_dimension : int, default 0
    filtration : np.ndarray, optional
    top_n : int, default 20

    Returns
    -------
    dict with keys:
        'features_disease', 'features_control' : pd.DataFrame
        'summary_disease',  'summary_control'  : pd.DataFrame
        'freq_disease',     'freq_control'     : pd.DataFrame
        'disease_specific', 'control_specific', 'shared' : list of str
        'comparison_df'    : pd.DataFrame — side-by-side gene comparison

    Example
    -------
        pph = GenPathHomology()
        betti_d, edges_d = pph.persistent_path_homology_from_digraph(
            X_disease, all_edges, 0, filtration
        )
        betti_c, edges_c = pph.persistent_path_homology_from_digraph(
            X_control, all_edges, 0, filtration
        )
        results = compare_group_features(
            betti_d, edges_d, betti_c, edges_c,
            n_nodes=len(gene_names), gene_names=gene_names,
        )
        print("Disease-specific genes:", results['disease_specific'])
        print(results['comparison_df'])
    """
    feat_d, summ_d = trace_persistent_features(
        betti_disease, edges_disease, n_nodes, gene_names,
        pph_instance, target_dimension, filtration
    )
    feat_c, summ_c = trace_persistent_features(
        betti_control, edges_control, n_nodes, gene_names,
        pph_instance, target_dimension, filtration
    )

    freq_d = summarise_gene_frequency(feat_d, top_n)
    freq_c = summarise_gene_frequency(feat_c, top_n)

    genes_d = set(freq_d["gene"]) if not freq_d.empty else set()
    genes_c = set(freq_c["gene"]) if not freq_c.empty else set()

    disease_specific = sorted(genes_d - genes_c)
    control_specific = sorted(genes_c - genes_d)
    shared           = sorted(genes_d & genes_c)

    all_genes   = sorted(genes_d | genes_c)
    freq_d_dict = dict(zip(freq_d["gene"], freq_d["n_steps"])) \
                  if not freq_d.empty else {}
    freq_c_dict = dict(zip(freq_c["gene"], freq_c["n_steps"])) \
                  if not freq_c.empty else {}

    comparison_records = []
    for gene in all_genes:
        nd = freq_d_dict.get(gene, 0)
        nc = freq_c_dict.get(gene, 0)
        comparison_records.append({
            "gene":            gene,
            "n_steps_disease": nd,
            "n_steps_control": nc,
            "difference":      nd - nc,   # positive = more active in disease
            "group": (
                "disease_specific" if nd > 0 and nc == 0 else
                "control_specific" if nc > 0 and nd == 0 else
                "shared"
            ),
        })

    comparison_df = (
        pd.DataFrame(comparison_records)
        .sort_values("difference", ascending=False)
        .reset_index(drop=True)
    )

    return {
        "features_disease":  feat_d,
        "features_control":  feat_c,
        "summary_disease":   summ_d,
        "summary_control":   summ_c,
        "freq_disease":      freq_d,
        "freq_control":      freq_c,
        "disease_specific":  disease_specific,
        "control_specific":  control_specific,
        "shared":            shared,
        "comparison_df":     comparison_df,
    }


# ══════════════════════════════════════════════════════════════════════════════
# Phase-level feature tracing
# ══════════════════════════════════════════════════════════════════════════════

def trace_features_by_phase(
    betti_series: list,
    edges_at_each_filtration: list,
    n_nodes: int,
    gene_names: List[str],
    pph_instance=None,
    target_dimension: int = 0,
    filtration: Optional[np.ndarray] = None,
    phases: Optional[Dict[str, Tuple[float, float]]] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Trace persistent features within each filtration phase.

    Splits the filtration into low/mid/high phases and identifies
    which genes drive topological features in each phase. Directly
    supports phase-based analysis in Paper 1 Pillar 2A:
        - Low phase (0.0–0.33): backbone genes at high correlation
        - Mid phase (0.33–0.66): bulk pathway structure genes
        - High phase (0.66–1.0): weakly correlated, likely noise

    Parameters
    ----------
    betti_series, edges_at_each_filtration, n_nodes, gene_names,
    pph_instance, target_dimension, filtration : see trace_persistent_features
    phases : dict, optional
        {phase_name: (low_threshold, high_threshold)}
        Default: {'low': (0.0, 0.33), 'mid': (0.33, 0.66),
                  'high': (0.66, 1.0)}

    Returns
    -------
    dict mapping phase_name -> pd.DataFrame of features in that phase

    Example
    -------
        phase_results = trace_features_by_phase(
            betti_0, edges_filt, len(gene_names), gene_names,
            filtration=filtration, target_dimension=0
        )
        print("Backbone genes (low phase):")
        print(phase_results['low'][['filtration_value','gene_members']])
        print("Bulk genes (mid phase):")
        print(phase_results['mid'][['filtration_value','gene_members']])
    """
    if phases is None:
        phases = {
            "low":  (0.00, 0.33),
            "mid":  (0.33, 0.66),
            "high": (0.66, 1.00),
        }

    if filtration is None:
        filtration = np.linspace(0, 1, len(betti_series))

    features_df, _ = trace_persistent_features(
        betti_series, edges_at_each_filtration,
        n_nodes, gene_names, pph_instance,
        target_dimension, filtration
    )

    phase_results = {}
    for name, (lo, hi) in phases.items():
        if features_df.empty:
            phase_results[name] = pd.DataFrame()
        else:
            mask = (
                (features_df["filtration_value"] >= lo) &
                (features_df["filtration_value"] <= hi)
            )
            phase_results[name] = features_df[mask].reset_index(drop=True)

    return phase_results
