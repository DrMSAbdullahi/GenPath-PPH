"""
genpath_pph/api.py
==================
High-level one-liner API for GenPath-PPH.

This module wraps the existing GenPathHomology, PathwayDataProcessor,
and statistical functions into a single clean interface so that a full
analysis can be run in as few lines as possible — similar to how other python packages work.

All computation is done by the existing core classes. This file only
chains them together and returns tidy results.

Quick start
-----------
    import numpy as np
    from genpath_pph.api import run_pathway, run_batch, GenPathAnalysis

    # ── Single pathway, one call ──────────────────────────────────────────
    # expression : np.ndarray (n_genes, n_samples)  — already pathway-subset
    # adj        : np.ndarray (n_genes, n_genes)    — directed adjacency

    result = run_pathway(X_disease, X_control, adj)
    print(result)
    # PathwayResult(significant=True, ks_pval_b0=0.012, cohend_b0=0.81, ...)

    # ── Batch over many pathways ──────────────────────────────────────────
    results_df = run_batch(
        pathway_ids      = select_path_ids,
        adj_matrices     = adj_matrices,          # dict: id -> DataFrame
        pathway_exprs    = pathways_expressions,  # dict: id -> DataFrame
        class_size       = 17,                    # samples per group
        n_permutations   = 5000,
    )
    # returns a tidy pandas DataFrame, one row per pathway

    # ── Object-oriented style ─────────────────────────────────────────────
    model = GenPathAnalysis(n_steps=100, n_permutations=1000)
    model.fit(X_disease, X_control, adj)

    b0_disease, b1_disease = model.betti_series("disease")
    b0_control, b1_control = model.betti_series("control")
    delta0, delta1         = model.delta_betti()
    phases                 = model.phase_summary()
    result                 = model.test()
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ── import from your existing modules ────────────────────────────────────────
from .core import GenPathHomology
from .utils import perform_ks_and_effectsize_tests


# ══════════════════════════════════════════════════════════════════════════════
# Result dataclass
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class PathwayResult:
    """
    Result of a single-pathway GenPath-PPH analysis.

    Attributes
    ----------
    ks_stat_b0, ks_stat_b1     : KS statistic for β₀ and β₁ series.
    cohend_b0, cohend_b1       : Cohen's d effect size.
    perm_pval_b0, perm_pval_b1 : Empirical p-value from permutation test.
    mean_diff_b0, mean_diff_b1 : Mean Betti difference (disease − control).
    significant                : True if BOTH dimensions pass KS + Cohen d.
    """
    ks_stat_b0:   float
    perm_pval_b0: float
    cohend_b0:    float
    mean_diff_b0: float
    ks_stat_b1:   float
    perm_pval_b1: float
    cohend_b1:    float
    mean_diff_b1: float
    significant:  bool
    alpha:        float = 0.05
    cohend_threshold: float = 0.5

    def __repr__(self):
        sig = "SIGNIFICANT" if self.significant else "not significant"
        return (
            f"PathwayResult({sig})\n"
            f"  β₀ → KS={self.ks_stat_b0:.4f}  perm-p={self.perm_pval_b0:.4f}"
            f"  Cohen-d={self.cohend_b0:.4f}  mean-diff={self.mean_diff_b0:.4f}\n"
            f"  β₁ → KS={self.ks_stat_b1:.4f}  perm-p={self.perm_pval_b1:.4f}"
            f"  Cohen-d={self.cohend_b1:.4f}  mean-diff={self.mean_diff_b1:.4f}"
        )

    def to_dict(self) -> dict:
        return {
            "ks_stat_b0":   self.ks_stat_b0,
            "perm_pval_b0": self.perm_pval_b0,
            "cohend_b0":    self.cohend_b0,
            "mean_diff_b0": self.mean_diff_b0,
            "ks_stat_b1":   self.ks_stat_b1,
            "perm_pval_b1": self.perm_pval_b1,
            "cohend_b1":    self.cohend_b1,
            "mean_diff_b1": self.mean_diff_b1,
            "significant":  self.significant,
        }


# ══════════════════════════════════════════════════════════════════════════════
# Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

def _adj_df_to_edges(adj_df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert an adjacency DataFrame (output of PathwayDataProcessor.get_adjacency_matrix)
    to a (rows, cols) edge index array compatible with persistent_path_homology_from_digraph.

    Parameters
    ----------
    adj_df : pd.DataFrame — square binary adjacency matrix with gene labels as index/columns.

    Returns
    -------
    edges : np.ndarray, shape (n_edges, 2) — integer row/col indices of edges.
    genes : list of str — gene labels in row/column order.
    """
    genes = list(adj_df.index)
    adj_np = adj_df.values.astype(int)
    rows, cols = np.nonzero(adj_np)
    edges = np.column_stack([rows, cols])
    return edges, genes


def _expr_df_to_array(expr_df: pd.DataFrame, genes: list) -> np.ndarray:
    """
    Reindex expression DataFrame to match gene order in adjacency matrix.

    Parameters
    ----------
    expr_df : pd.DataFrame — gene expression (genes × samples).
    genes   : list — gene order from adjacency matrix.

    Returns
    -------
    np.ndarray, shape (n_genes, n_samples)
    """
    common = [g for g in genes if g in expr_df.index]
    return expr_df.loc[common].values


def _run_pph_single(
    X: np.ndarray,
    edges: np.ndarray,
    target_dim: int,
    filtration: np.ndarray,
    distance_type: str,
) -> list:
    """Run PPH for one group and one dimension."""
    pph = GenPathHomology()
    betti, _ = pph.persistent_path_homology_from_digraph(
        cloudpoints=X,
        all_edges=edges,
        target_dimension=target_dim,
        filtration=filtration,
        distance_type=distance_type,
    )
    return betti


def _stat_single_dim(
    betti_control: list,
    betti_disease: list,
    n_permutations: int,
    seed: int,
) -> Tuple[float, float, float]:
    """
    Run KS + Cohen d permutation tests for one Betti dimension.
    Uses the existing perform_ks_and_effectsize_tests under the hood
    by packaging as a one-pathway dict.

    Returns (ks_stat, perm_pval, cohend, mean_diff)
    """
    # Reuse the existing permutation_test_effect_size directly
    from .utils import permutation_test_effect_size

    b_c = np.array(betti_control, dtype=float)
    b_d = np.array(betti_disease, dtype=float)

    ks_stat, _, ks_pval    = permutation_test_effect_size(b_c, b_d, n_permutations, "ks",      seed)
    cohend,  mean_d, cd_pval = permutation_test_effect_size(b_c, b_d, n_permutations, "cohen_d", seed)

    return float(ks_stat), float(ks_pval), float(cohend), float(mean_d) if not np.isnan(mean_d) else 0.0


# ══════════════════════════════════════════════════════════════════════════════
# One-liner function: single pathway
# ══════════════════════════════════════════════════════════════════════════════

def run_pathway(
    X_disease: np.ndarray,
    X_control: np.ndarray,
    adj: np.ndarray,
    filtration_scale: float = 1.0,
    step_size: float = 0.01,
    n_permutations: int = 1000,
    alpha: float = 0.05,
    cohend_threshold: float = 0.5,
    distance_type: str = "1-abs-correlation",
    seed: int = 0,
) -> PathwayResult:
    """
    Run the full GenPath-PPH pipeline on a single pathway in one call.

    Parameters
    ----------
    X_disease : np.ndarray, shape (n_genes, n_disease_samples)
        Gene expression for disease group, rows already subset to pathway genes
        and sorted to match rows/cols of adj.
    X_control : np.ndarray, shape (n_genes, n_control_samples)
        Gene expression for control group, same gene order as X_disease.
    adj : np.ndarray, shape (n_genes, n_genes)
        Binary directed adjacency matrix of the pathway graph.
        adj[i,j] = 1 means directed edge gene_i → gene_j.
        Can also be a pd.DataFrame (output of PathwayDataProcessor.get_adjacency_matrix).
    filtration_scale : float, default 1.0
        Maximum filtration value (distance range is [0, filtration_scale]).
    step_size : float, default 0.01
        Filtration step size (101 steps with scale=1, step=0.01).
    n_permutations : int, default 1000
        Number of permutations for empirical p-value.
    alpha : float, default 0.05
        Significance threshold.
    cohend_threshold : float, default 0.5
        Minimum |Cohen's d| for significance.
    distance_type : str, default '1-abs-correlation'
        Distance metric. Options: '1-abs-correlation', '1-correlation', 'euclidean'.
    seed : int, default 0
        Random seed for reproducibility.

    Returns
    -------
    PathwayResult

    Example
    -------
        import numpy as np
        from genpath_pph.api import run_pathway

        # X_disease, X_control: (n_genes, n_samples) expression arrays
        # adj: (n_genes, n_genes) binary adjacency matrix
        result = run_pathway(X_disease, X_control, adj)
        print(result)
        # PathwayResult(SIGNIFICANT)
        #   β₀ → KS=0.4231  perm-p=0.0020  Cohen-d=0.8813  mean-diff=2.3100
        #   β₁ → KS=0.3812  perm-p=0.0140  Cohen-d=0.7241  mean-diff=0.9800
    """
    # Accept both np.ndarray and pd.DataFrame adjacency
    if isinstance(adj, pd.DataFrame):
        adj_np = adj.values.astype(int)
    else:
        adj_np = np.asarray(adj, dtype=int)

    rows, cols = np.nonzero(adj_np)
    edges = np.column_stack([rows, cols]) if len(rows) > 0 else np.empty((0, 2), dtype=int)

    filtration = np.arange(0, filtration_scale, step_size)

    # β₀
    b0_c = _run_pph_single(X_control, edges, 0, filtration, distance_type)
    b0_d = _run_pph_single(X_disease, edges, 0, filtration, distance_type)

    # β₁
    b1_c = _run_pph_single(X_control, edges, 1, filtration, distance_type)
    b1_d = _run_pph_single(X_disease, edges, 1, filtration, distance_type)

    # Statistics
    ks0, pval0, cd0, md0 = _stat_single_dim(b0_c, b0_d, n_permutations, seed)
    ks1, pval1, cd1, md1 = _stat_single_dim(b1_c, b1_d, n_permutations, seed)

    sig0 = (pval0 < alpha) and (abs(cd0) >= cohend_threshold)
    sig1 = (pval1 < alpha) and (abs(cd1) >= cohend_threshold)

    return PathwayResult(
        ks_stat_b0=ks0,   perm_pval_b0=pval0, cohend_b0=cd0,   mean_diff_b0=md0,
        ks_stat_b1=ks1,   perm_pval_b1=pval1, cohend_b1=cd1,   mean_diff_b1=md1,
        significant=sig0 and sig1,
        alpha=alpha, cohend_threshold=cohend_threshold,
    )


# ══════════════════════════════════════════════════════════════════════════════
# One-liner function: batch over many pathways
# ══════════════════════════════════════════════════════════════════════════════

def run_batch(
    pathway_ids: List[str],
    adj_matrices: Dict[str, pd.DataFrame],
    pathway_exprs: Dict[str, pd.DataFrame],
    class_size: int = 17,
    filtration_scale: float = 1.0,
    step_size: float = 0.01,
    n_permutations: int = 5000,
    alpha: float = 0.05,
    cohend_threshold: float = 0.5,
    distance_type: str = "1-abs-correlation",
    seed: int = 0,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Run GenPath-PPH across multiple pathways in one call.

    This directly wraps GenPathHomology.compute_betti_numbers() for PPH
    computation and perform_ks_and_effectsize_tests() for statistics,
    giving you a single tidy DataFrame as output.

    Parameters
    ----------
    pathway_ids : list of str
        KEGG pathway IDs to analyse (e.g. ['hsa04151', 'hsa04115']).
    adj_matrices : dict
        Maps pathway_id -> adjacency DataFrame
        (output of PathwayDataProcessor.get_adjacency_matrix).
    pathway_exprs : dict
        Maps pathway_id -> expression DataFrame (genes × samples).
        Columns: first class_size samples = control, next class_size = disease.
        (output of extract_pathway_expressions in utils.py)
    class_size : int, default 17
        Number of samples per group (control and disease).
    filtration_scale : float, default 1.0
    step_size : float, default 0.01
    n_permutations : int, default 5000
    alpha : float, default 0.05
    cohend_threshold : float, default 0.5
    distance_type : str, default '1-abs-correlation'
    seed : int, default 0
    verbose : bool, default True

    Returns
    -------
    pd.DataFrame with columns:
        path_id, mean_diff_b0, es_b0, perm_pval_b0, fdr_pval_b0,
                 mean_diff_b1, es_b1, perm_pval_b1, fdr_pval_b1,
        significant (True if both dimensions pass after FDR correction)

    Example
    -------
        from genpath_pph.api import run_batch
        from genpath_pph.utils import extract_adjacency_matrices, extract_pathway_expressions

        adj_matrices    = extract_adjacency_matrices(select_path_ids, ...)
        pathway_exprs   = extract_pathway_expressions(select_path_ids, ...)

        results = run_batch(
            pathway_ids   = select_path_ids,
            adj_matrices  = adj_matrices,
            pathway_exprs = pathway_exprs,
            class_size    = 17,
        )
        significant = results[results["significant"]]
        print(f"{len(significant)} significant pathways found")
    """
    pph = GenPathHomology()

    # ── Step 1: compute Betti series for all pathways (both dims) ─────────
    if verbose:
        print(f"Computing Betti series for {len(pathway_ids)} pathways...")

    betti_dim0 = pph.compute_betti_numbers(
        select_path_ids    = pathway_ids,
        adj_matrices       = {k: v.values.astype(int) for k, v in adj_matrices.items()},
        pathways_expressions = pathway_exprs,
        target_dimension   = 0,
        filtration_scale   = filtration_scale,
        step_size          = step_size,
        class_size         = class_size,
        distance_type      = distance_type,
    )

    betti_dim1 = pph.compute_betti_numbers(
        select_path_ids    = pathway_ids,
        adj_matrices       = {k: v.values.astype(int) for k, v in adj_matrices.items()},
        pathways_expressions = pathway_exprs,
        target_dimension   = 1,
        filtration_scale   = filtration_scale,
        step_size          = step_size,
        class_size         = class_size,
        distance_type      = distance_type,
    )

    # ── Step 2: repackage into the dict format expected by perform_ks_and_effectsize_tests ──
    # That function expects: betti_dict[pid]['control'] and betti_dict[pid]['disease']
    betti_dict_b0 = {
        pid: {
            "control": np.array(betti_dim0["control"][pid], dtype=float),
            "disease": np.array(betti_dim0["disease"][pid], dtype=float),
        }
        for pid in pathway_ids
    }
    betti_dict_b1 = {
        pid: {
            "control": np.array(betti_dim1["control"][pid], dtype=float),
            "disease": np.array(betti_dim1["disease"][pid], dtype=float),
        }
        for pid in pathway_ids
    }

    # ── Step 3: statistical tests using the existing function ─────────────
    if verbose:
        print("Running statistical tests (β₀)...")
    df_b0 = perform_ks_and_effectsize_tests(
        betti_dict     = betti_dict_b0,
        path_ids       = pathway_ids,
        num_permutations = n_permutations,
        seed           = seed,
    ).rename(columns={
        "mean_diff_observed": "mean_diff_b0",
        "es_observed":        "es_b0",
        "es_raw_pvalue":      "es_raw_pval_b0",
        "ks_observed":        "ks_stat_b0",
        "ks_raw_pvalue":      "ks_raw_pval_b0",
        "es_corrected_pvalue":"fdr_es_b0",
        "ks_corrected_pvalue":"fdr_ks_b0",
    })

    if verbose:
        print("Running statistical tests (β₁)...")
    df_b1 = perform_ks_and_effectsize_tests(
        betti_dict     = betti_dict_b1,
        path_ids       = pathway_ids,
        num_permutations = n_permutations,
        seed           = seed,
    ).rename(columns={
        "mean_diff_observed": "mean_diff_b1",
        "es_observed":        "es_b1",
        "es_raw_pvalue":      "es_raw_pval_b1",
        "ks_observed":        "ks_stat_b1",
        "ks_raw_pvalue":      "ks_raw_pval_b1",
        "es_corrected_pvalue":"fdr_es_b1",
        "ks_corrected_pvalue":"fdr_ks_b1",
    })

    # ── Step 4: merge and add significance flag ────────────────────────────
    results = df_b0.merge(df_b1, on="path_id")

    results["significant"] = (
        (results["fdr_ks_b0"] < alpha) & (results["es_b0"].abs() >= cohend_threshold) &
        (results["fdr_ks_b1"] < alpha) & (results["es_b1"].abs() >= cohend_threshold)
    )

    if verbose:
        n_sig = results["significant"].sum()
        print(f"\nDone. {n_sig}/{len(pathway_ids)} pathways significant.")

    return results


# ══════════════════════════════════════════════════════════════════════════════
# Object-oriented interface
# ══════════════════════════════════════════════════════════════════════════════

class GenPathAnalysis:
    """
    Object-oriented interface to GenPath-PPH for a single pathway.

    Mirrors the scikit-learn fit/predict pattern. Internally delegates
    all computation to GenPathHomology.

    Parameters
    ----------
    filtration_scale : float, default 1.0
    step_size : float, default 0.01
    n_permutations : int, default 1000
    alpha : float, default 0.05
    cohend_threshold : float, default 0.5
    distance_type : str, default '1-abs-correlation'
    seed : int, default 0

    Example
    -------
        model = GenPathAnalysis(n_permutations=1000)
        model.fit(X_disease, X_control, adj)

        b0_d, b1_d = model.betti_series("disease")
        b0_c, b1_c = model.betti_series("control")
        delta0, delta1 = model.delta_betti()
        phases = model.phase_summary()
        result = model.test()
        print(result)
    """

    def __init__(
        self,
        filtration_scale: float = 1.0,
        step_size: float = 0.01,
        n_permutations: int = 1000,
        alpha: float = 0.05,
        cohend_threshold: float = 0.5,
        distance_type: str = "1-abs-correlation",
        seed: int = 0,
    ):
        self.filtration_scale  = filtration_scale
        self.step_size         = step_size
        self.n_permutations    = n_permutations
        self.alpha             = alpha
        self.cohend_threshold  = cohend_threshold
        self.distance_type     = distance_type
        self.seed              = seed

        self._b0_disease: Optional[np.ndarray] = None
        self._b1_disease: Optional[np.ndarray] = None
        self._b0_control: Optional[np.ndarray] = None
        self._b1_control: Optional[np.ndarray] = None
        self._thresholds: Optional[np.ndarray] = None
        self._fitted = False

    def fit(
        self,
        X_disease: np.ndarray,
        X_control: np.ndarray,
        adj,               # np.ndarray or pd.DataFrame
    ) -> "GenPathAnalysis":
        """
        Compute PPH Betti series for disease and control groups.

        Parameters
        ----------
        X_disease : np.ndarray, shape (n_genes, n_disease_samples)
        X_control : np.ndarray, shape (n_genes, n_control_samples)
        adj : np.ndarray or pd.DataFrame, shape (n_genes, n_genes)
            Binary directed adjacency matrix.

        Returns
        -------
        self
        """
        if isinstance(adj, pd.DataFrame):
            adj_np = adj.values.astype(int)
        else:
            adj_np = np.asarray(adj, dtype=int)

        rows, cols = np.nonzero(adj_np)
        edges = np.column_stack([rows, cols]) if len(rows) > 0 else np.empty((0,2), dtype=int)

        self._thresholds = np.arange(0, self.filtration_scale, self.step_size)

        self._b0_control = np.array(_run_pph_single(X_control, edges, 0, self._thresholds, self.distance_type))
        self._b0_disease = np.array(_run_pph_single(X_disease, edges, 0, self._thresholds, self.distance_type))
        self._b1_control = np.array(_run_pph_single(X_control, edges, 1, self._thresholds, self.distance_type))
        self._b1_disease = np.array(_run_pph_single(X_disease, edges, 1, self._thresholds, self.distance_type))

        self._fitted = True
        return self

    def betti_series(self, group: str = "disease") -> Tuple[np.ndarray, np.ndarray]:
        """
        Return (β₀, β₁) Betti series for the requested group.

        Parameters
        ----------
        group : str — 'disease' or 'control'

        Returns
        -------
        beta0 : np.ndarray, shape (n_steps,)
        beta1 : np.ndarray, shape (n_steps,)
        """
        self._check_fitted()
        if group == "disease":
            return self._b0_disease.copy(), self._b1_disease.copy()
        elif group == "control":
            return self._b0_control.copy(), self._b1_control.copy()
        else:
            raise ValueError(f"group must be 'disease' or 'control', got '{group}'.")

    def delta_betti(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Return difference curves Δβ₀(t) and Δβ₁(t) = disease − control.

        Returns
        -------
        delta0 : np.ndarray
        delta1 : np.ndarray
        """
        self._check_fitted()
        return self._b0_disease - self._b0_control, self._b1_disease - self._b1_control

    def thresholds(self) -> np.ndarray:
        """Return the filtration threshold values."""
        self._check_fitted()
        return self._thresholds.copy()

    def aggregate_score(self) -> Dict[str, float]:
        """Mean Δβ across filtration (pathway dysregulation score)."""
        d0, d1 = self.delta_betti()
        return {"delta_beta0": float(d0.mean()), "delta_beta1": float(d1.mean())}

    def phase_summary(
        self,
        phases: Optional[Dict[str, Tuple[float, float]]] = None,
    ) -> Dict[str, Dict[str, float]]:
        """
        Mean Δβ within each filtration phase.

        Parameters
        ----------
        phases : dict, optional
            {phase_name: (low_threshold, high_threshold)}
            Default: {'low': (0.0, 0.33), 'mid': (0.33, 0.66), 'high': (0.66, 1.0)}

        Returns
        -------
        dict: {phase_name: {'delta_beta0': float, 'delta_beta1': float}}
        """
        self._check_fitted()
        if phases is None:
            phases = {"low": (0.0, 0.33), "mid": (0.33, 0.66), "high": (0.66, 1.0)}
        d0, d1 = self.delta_betti()
        summary = {}
        for name, (lo, hi) in phases.items():
            mask = (self._thresholds >= lo) & (self._thresholds <= hi)
            summary[name] = {
                "delta_beta0": float(d0[mask].mean()) if mask.any() else 0.0,
                "delta_beta1": float(d1[mask].mean()) if mask.any() else 0.0,
            }
        return summary

    def test(self) -> PathwayResult:
        """
        Run KS + Cohen d + permutation test on the fitted Betti series.

        Returns
        -------
        PathwayResult
        """
        self._check_fitted()

        ks0, pval0, cd0, md0 = _stat_single_dim(
            self._b0_control.tolist(), self._b0_disease.tolist(),
            self.n_permutations, self.seed
        )
        ks1, pval1, cd1, md1 = _stat_single_dim(
            self._b1_control.tolist(), self._b1_disease.tolist(),
            self.n_permutations, self.seed
        )

        sig0 = (pval0 < self.alpha) and (abs(cd0) >= self.cohend_threshold)
        sig1 = (pval1 < self.alpha) and (abs(cd1) >= self.cohend_threshold)

        return PathwayResult(
            ks_stat_b0=ks0,   perm_pval_b0=pval0, cohend_b0=cd0,   mean_diff_b0=md0,
            ks_stat_b1=ks1,   perm_pval_b1=pval1, cohend_b1=cd1,   mean_diff_b1=md1,
            significant=sig0 and sig1,
            alpha=self.alpha, cohend_threshold=self.cohend_threshold,
        )

    def _check_fitted(self):
        if not self._fitted:
            raise RuntimeError("Call fit() before accessing results.")
