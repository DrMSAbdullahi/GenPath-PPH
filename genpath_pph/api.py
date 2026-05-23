"""
genpath_pph/api.py
==================
High-level one-liner API for GenPath-PPH.

This module wraps the existing GenPathHomology, PathwayDataProcessor,
and statistical functions into a single clean interface so that a full
analysis can be run in as few lines as possible — similar to how other
Python packages such as scikit-learn or GUDHI work.

All computation is done by the existing core classes. This file only
chains them together and returns tidy results.

Significance criteria (consistent with Abdullahi et al., 2025, CSBJ)
---------------------------------------------------------------------
A pathway is significant if ALL FOUR corrected p-values < alpha:
    - KS test permutation p-value for β₀
    - KS test permutation p-value for β₁
    - Cohen's d permutation p-value for β₀
    - Cohen's d permutation p-value for β₁

For run_pathway() (single pathway), raw permutation p-values are used
since there is only one pathway and no multiple testing correction is
needed. For run_batch() (many pathways), Benjamini-Hochberg FDR
correction is applied across all pathways before significance is
determined, exactly as in the CSBJ paper.

Quick start
-----------
    import numpy as np
    from genpath_pph.api import run_pathway, run_batch, GenPathAnalysis

    # ── Single pathway, one call ──────────────────────────────────────────
    result = run_pathway(X_disease, X_control, adj)
    print(result)
    # PathwayResult(SIGNIFICANT)
    #   β₀ → KS=0.20  KS-perm-p=0.000  Cohen-d=-0.64  Cohen-perm-p=0.021  mean-diff=-0.40
    #   β₁ → KS=0.72  KS-perm-p=0.000  Cohen-d=2.10   Cohen-perm-p=0.000  mean-diff=0.72

    # ── Batch over many pathways ──────────────────────────────────────────
    results_df = run_batch(
        pathway_ids   = select_path_ids,
        adj_matrices  = adj_matrices,
        pathway_exprs = pathways_expressions,
        class_size    = 17,
        n_permutations = 5000,
    )

    # ── Object-oriented style ─────────────────────────────────────────────
    model = GenPathAnalysis(n_permutations=1000)
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

# ── import from existing modules ─────────────────────────────────────────────
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
    ks_stat_b0, ks_stat_b1         : KS statistic for β₀ and β₁ series.
    perm_pval_b0, perm_pval_b1     : Permutation p-value for KS test.
    cohend_b0, cohend_b1           : Cohen's d effect size.
    cohend_pval_b0, cohend_pval_b1 : Permutation p-value for Cohen's d.
    mean_diff_b0, mean_diff_b1     : Mean Betti difference (control − disease).
    significant                    : True if ALL FOUR p-values < alpha.
    alpha                          : Significance threshold used.

    Notes
    -----
    Significance follows Abdullahi et al. (2025):
        significant if KS-perm-p < alpha AND Cohen-perm-p < alpha
        for BOTH β₀ and β₁ (four conditions total).
    mean_diff is defined as mean(control) − mean(disease) across filtration.
    """
    ks_stat_b0:     float
    perm_pval_b0:   float
    cohend_b0:      float
    cohend_pval_b0: float
    mean_diff_b0:   float
    ks_stat_b1:     float
    perm_pval_b1:   float
    cohend_b1:      float
    cohend_pval_b1: float
    mean_diff_b1:   float
    significant:    bool
    alpha:          float = 0.05

    def __repr__(self):
        sig = "SIGNIFICANT" if self.significant else "not significant"
        return (
            f"PathwayResult({sig})\n"
            f"  β₀ → KS={self.ks_stat_b0:.4f}  KS-perm-p={self.perm_pval_b0:.4f}"
            f"  Cohen-d={self.cohend_b0:.4f}  Cohen-perm-p={self.cohend_pval_b0:.4f}"
            f"  mean-diff={self.mean_diff_b0:.4f}\n"
            f"  β₁ → KS={self.ks_stat_b1:.4f}  KS-perm-p={self.perm_pval_b1:.4f}"
            f"  Cohen-d={self.cohend_b1:.4f}  Cohen-perm-p={self.cohend_pval_b1:.4f}"
            f"  mean-diff={self.mean_diff_b1:.4f}"
        )

    def to_dict(self) -> dict:
        return {
            "ks_stat_b0":     self.ks_stat_b0,
            "perm_pval_b0":   self.perm_pval_b0,
            "cohend_b0":      self.cohend_b0,
            "cohend_pval_b0": self.cohend_pval_b0,
            "mean_diff_b0":   self.mean_diff_b0,
            "ks_stat_b1":     self.ks_stat_b1,
            "perm_pval_b1":   self.perm_pval_b1,
            "cohend_b1":      self.cohend_b1,
            "cohend_pval_b1": self.cohend_pval_b1,
            "mean_diff_b1":   self.mean_diff_b1,
            "significant":    self.significant,
        }


# ══════════════════════════════════════════════════════════════════════════════
# Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

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
) -> Tuple[float, float, float, float, float]:
    """
    Run KS + Cohen d permutation tests for one Betti dimension.

    Returns
    -------
    ks_stat, ks_pval, cohend, cohend_pval, mean_diff
    """
    from .utils import permutation_test_effect_size

    b_c = np.array(betti_control, dtype=float)
    b_d = np.array(betti_disease, dtype=float)

    # KS permutation test
    ks_stat, _, ks_pval = permutation_test_effect_size(
        b_c, b_d, n_permutations, "ks", seed
    )

    # Cohen's d permutation test
    cohend, mean_d, cd_pval = permutation_test_effect_size(
        b_c, b_d, n_permutations, "cohen_d", seed
    )

    mean_diff = float(mean_d) if not np.isnan(mean_d) else 0.0

    return (
        float(ks_stat),
        float(ks_pval),
        float(cohend),
        float(cd_pval),
        mean_diff,
    )


# ══════════════════════════════════════════════════════════════════════════════
# One-liner function: single pathway
# ══════════════════════════════════════════════════════════════════════════════

def run_pathway(
    X_disease: np.ndarray,
    X_control: np.ndarray,
    adj,
    filtration_scale: float = 1.0,
    step_size: float = 0.01,
    n_permutations: int = 1000,
    alpha: float = 0.05,
    distance_type: str = "1-abs-correlation",
    seed: int = 0,
) -> PathwayResult:
    """
    Run the full GenPath-PPH pipeline on a single pathway in one call.

    Parameters
    ----------
    X_disease : np.ndarray, shape (n_genes, n_disease_samples)
        Gene expression for disease group, rows subset to pathway genes.
    X_control : np.ndarray, shape (n_genes, n_control_samples)
        Gene expression for control group, same gene order as X_disease.
    adj : np.ndarray or pd.DataFrame, shape (n_genes, n_genes)
        Binary directed adjacency matrix of the pathway graph.
        adj[i,j] = 1 means directed edge gene_i → gene_j.
    filtration_scale : float, default 1.0
        Maximum filtration value.
    step_size : float, default 0.01
        Filtration step size.
    n_permutations : int, default 1000
        Number of permutations for empirical p-values.
    alpha : float, default 0.05
        Significance threshold for all four p-values.
    distance_type : str, default '1-abs-correlation'
        Distance metric. Options: '1-abs-correlation', '1-correlation',
        'euclidean'.
    seed : int, default 0
        Random seed for reproducibility.

    Returns
    -------
    PathwayResult

    Notes
    -----
    Significance follows Abdullahi et al. (2025, CSBJ):
    A pathway is significant if ALL FOUR p-values < alpha:
        KS-perm-p (β₀), KS-perm-p (β₁),
        Cohen-perm-p (β₀), Cohen-perm-p (β₁).
    Raw permutation p-values are used here (single pathway, no multiple
    testing correction needed). For multi-pathway analysis use run_batch()
    which applies Benjamini-Hochberg FDR correction.

    Example
    -------
        from genpath_pph import run_pathway
        result = run_pathway(X_disease, X_control, adj_matrix)
        print(result)
    """
    # Accept both np.ndarray and pd.DataFrame adjacency
    if isinstance(adj, pd.DataFrame):
        adj_np = adj.values.astype(int)
    else:
        adj_np = np.asarray(adj, dtype=int)

    rows, cols = np.nonzero(adj_np)
    edges = (np.column_stack([rows, cols]) if len(rows) > 0
             else np.empty((0, 2), dtype=int))

    filtration = np.arange(0, filtration_scale, step_size)

    # Compute Betti series
    b0_c = _run_pph_single(X_control, edges, 0, filtration, distance_type)
    b0_d = _run_pph_single(X_disease, edges, 0, filtration, distance_type)
    b1_c = _run_pph_single(X_control, edges, 1, filtration, distance_type)
    b1_d = _run_pph_single(X_disease, edges, 1, filtration, distance_type)

    # Statistics — note: control vs disease order matches mean_diff definition
    ks0, ksp0, cd0, cdp0, md0 = _stat_single_dim(b0_c, b0_d, n_permutations, seed)
    ks1, ksp1, cd1, cdp1, md1 = _stat_single_dim(b1_c, b1_d, n_permutations, seed)

    # Significance: all four p-values must be < alpha
    sig0 = (ksp0 < alpha) and (cdp0 < alpha)
    sig1 = (ksp1 < alpha) and (cdp1 < alpha)

    return PathwayResult(
        ks_stat_b0=ks0,   perm_pval_b0=ksp0,
        cohend_b0=cd0,    cohend_pval_b0=cdp0,  mean_diff_b0=md0,
        ks_stat_b1=ks1,   perm_pval_b1=ksp1,
        cohend_b1=cd1,    cohend_pval_b1=cdp1,  mean_diff_b1=md1,
        significant=sig0 and sig1,
        alpha=alpha,
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
    distance_type: str = "1-abs-correlation",
    seed: int = 0,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Run GenPath-PPH across multiple pathways in one call.

    Wraps GenPathHomology.compute_betti_numbers() for PPH computation
    and perform_ks_and_effectsize_tests() for statistics, applying
    Benjamini-Hochberg FDR correction exactly as in Abdullahi et al.
    (2025, CSBJ). A pathway is significant if all four FDR-corrected
    p-values (KS and Cohen's d for both β₀ and β₁) are < alpha.

    Parameters
    ----------
    pathway_ids : list of str
        KEGG pathway IDs to analyse (e.g. ['hsa04151', 'hsa04115']).
    adj_matrices : dict
        Maps pathway_id -> adjacency DataFrame
        (output of PathwayDataProcessor.get_adjacency_matrix).
    pathway_exprs : dict
        Maps pathway_id -> expression DataFrame (genes x samples).
        Columns: first class_size columns = control,
                 next class_size columns = disease.
        (output of extract_pathway_expressions in utils.py)
    class_size : int, default 17
        Number of samples per group.
    filtration_scale : float, default 1.0
    step_size : float, default 0.01
    n_permutations : int, default 5000
    alpha : float, default 0.05
    distance_type : str, default '1-abs-correlation'
    seed : int, default 0
    verbose : bool, default True

    Returns
    -------
    pd.DataFrame with columns:
        path_id,
        mean_diff_b0, es_b0, es_raw_pval_b0, fdr_es_b0,
                      ks_stat_b0, ks_raw_pval_b0, fdr_ks_b0,
        mean_diff_b1, es_b1, es_raw_pval_b1, fdr_es_b1,
                      ks_stat_b1, ks_raw_pval_b1, fdr_ks_b1,
        significant

    Example
    -------
        from genpath_pph import run_batch
        from genpath_pph import extract_adjacency_matrices
        from genpath_pph import extract_pathway_expressions

        adj_matrices  = extract_adjacency_matrices(select_path_ids, ...)
        pathway_exprs = extract_pathway_expressions(select_path_ids, ...)

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

    # ── Step 1: compute Betti series for all pathways ─────────────────────
    if verbose:
        print(f"Computing Betti series for {len(pathway_ids)} pathways...")

    betti_dim0 = pph.compute_betti_numbers(
        select_path_ids      = pathway_ids,
        adj_matrices         = {k: v.values.astype(int)
                                for k, v in adj_matrices.items()},
        pathways_expressions = pathway_exprs,
        target_dimension     = 0,
        filtration_scale     = filtration_scale,
        step_size            = step_size,
        class_size           = class_size,
        distance_type        = distance_type,
    )

    betti_dim1 = pph.compute_betti_numbers(
        select_path_ids      = pathway_ids,
        adj_matrices         = {k: v.values.astype(int)
                                for k, v in adj_matrices.items()},
        pathways_expressions = pathway_exprs,
        target_dimension     = 1,
        filtration_scale     = filtration_scale,
        step_size            = step_size,
        class_size           = class_size,
        distance_type        = distance_type,
    )

    # ── Step 2: repackage into dict format for statistical tests ──────────
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

    # ── Step 3: statistical tests with BH correction ──────────────────────
    if verbose:
        print("Running statistical tests (β₀)...")
    df_b0 = perform_ks_and_effectsize_tests(
        betti_dict       = betti_dict_b0,
        path_ids         = pathway_ids,
        num_permutations = n_permutations,
        seed             = seed,
    ).rename(columns={
        "mean_diff_observed":  "mean_diff_b0",
        "es_observed":         "es_b0",
        "es_raw_pvalue":       "es_raw_pval_b0",
        "ks_observed":         "ks_stat_b0",
        "ks_raw_pvalue":       "ks_raw_pval_b0",
        "es_corrected_pvalue": "fdr_es_b0",
        "ks_corrected_pvalue": "fdr_ks_b0",
    })

    if verbose:
        print("Running statistical tests (β₁)...")
    df_b1 = perform_ks_and_effectsize_tests(
        betti_dict       = betti_dict_b1,
        path_ids         = pathway_ids,
        num_permutations = n_permutations,
        seed             = seed,
    ).rename(columns={
        "mean_diff_observed":  "mean_diff_b1",
        "es_observed":         "es_b1",
        "es_raw_pvalue":       "es_raw_pval_b1",
        "ks_observed":         "ks_stat_b1",
        "ks_raw_pvalue":       "ks_raw_pval_b1",
        "es_corrected_pvalue": "fdr_es_b1",
        "ks_corrected_pvalue": "fdr_ks_b1",
    })

    # ── Step 4: merge and add significance flag ────────────────────────────
    results = df_b0.merge(df_b1, on="path_id")

    # Significance: all four FDR-corrected p-values < alpha
    # (KS and Cohen's d for both β₀ and β₁)
    results["significant"] = (
        (results["fdr_ks_b0"] < alpha) &
        (results["fdr_es_b0"] < alpha) &
        (results["fdr_ks_b1"] < alpha) &
        (results["fdr_es_b1"] < alpha)
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

    Significance follows Abdullahi et al. (2025, CSBJ): a pathway is
    significant if ALL FOUR permutation p-values < alpha (KS and
    Cohen's d for both β₀ and β₁).

    Parameters
    ----------
    filtration_scale : float, default 1.0
    step_size : float, default 0.01
    n_permutations : int, default 1000
    alpha : float, default 0.05
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
        distance_type: str = "1-abs-correlation",
        seed: int = 0,
    ):
        self.filtration_scale = filtration_scale
        self.step_size        = step_size
        self.n_permutations   = n_permutations
        self.alpha            = alpha
        self.distance_type    = distance_type
        self.seed             = seed

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
        adj,
    ) -> "GenPathAnalysis":
        """
        Compute PPH Betti series for disease and control groups.

        Parameters
        ----------
        X_disease : np.ndarray, shape (n_genes, n_disease_samples)
        X_control : np.ndarray, shape (n_genes, n_control_samples)
        adj : np.ndarray or pd.DataFrame, shape (n_genes, n_genes)

        Returns
        -------
        self
        """
        if isinstance(adj, pd.DataFrame):
            adj_np = adj.values.astype(int)
        else:
            adj_np = np.asarray(adj, dtype=int)

        rows, cols = np.nonzero(adj_np)
        edges = (np.column_stack([rows, cols]) if len(rows) > 0
                 else np.empty((0, 2), dtype=int))

        self._thresholds = np.arange(0, self.filtration_scale, self.step_size)

        self._b0_control = np.array(
            _run_pph_single(X_control, edges, 0, self._thresholds, self.distance_type))
        self._b0_disease = np.array(
            _run_pph_single(X_disease, edges, 0, self._thresholds, self.distance_type))
        self._b1_control = np.array(
            _run_pph_single(X_control, edges, 1, self._thresholds, self.distance_type))
        self._b1_disease = np.array(
            _run_pph_single(X_disease, edges, 1, self._thresholds, self.distance_type))

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
            raise ValueError(
                f"group must be 'disease' or 'control', got '{group}'."
            )

    def delta_betti(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Return difference curves Δβ₀(t) and Δβ₁(t) = control − disease.

        Note: control − disease, consistent with mean_diff definition
        in Abdullahi et al. (2025).

        Returns
        -------
        delta0 : np.ndarray
        delta1 : np.ndarray
        """
        self._check_fitted()
        return (self._b0_control - self._b0_disease,
                self._b1_control - self._b1_disease)

    def thresholds(self) -> np.ndarray:
        """Return the filtration threshold values."""
        self._check_fitted()
        return self._thresholds.copy()

    def aggregate_score(self) -> Dict[str, float]:
        """
        Mean Δβ across filtration (pathway dysregulation score).
        Defined as mean(control) − mean(disease).
        """
        d0, d1 = self.delta_betti()
        return {
            "delta_beta0": float(d0.mean()),
            "delta_beta1": float(d1.mean()),
        }

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
            Default: {'low': (0.0, 0.33),
                      'mid': (0.33, 0.66),
                      'high': (0.66, 1.0)}

        Returns
        -------
        dict: {phase_name: {'delta_beta0': float, 'delta_beta1': float}}
        """
        self._check_fitted()
        if phases is None:
            phases = {
                "low":  (0.00, 0.33),
                "mid":  (0.33, 0.66),
                "high": (0.66, 1.00),
            }
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
        Run KS + Cohen d permutation tests on the fitted Betti series.

        Significance follows Abdullahi et al. (2025, CSBJ):
        significant if ALL FOUR permutation p-values < alpha.

        Returns
        -------
        PathwayResult
        """
        self._check_fitted()

        ks0, ksp0, cd0, cdp0, md0 = _stat_single_dim(
            self._b0_control.tolist(),
            self._b0_disease.tolist(),
            self.n_permutations,
            self.seed,
        )
        ks1, ksp1, cd1, cdp1, md1 = _stat_single_dim(
            self._b1_control.tolist(),
            self._b1_disease.tolist(),
            self.n_permutations,
            self.seed,
        )

        sig0 = (ksp0 < self.alpha) and (cdp0 < self.alpha)
        sig1 = (ksp1 < self.alpha) and (cdp1 < self.alpha)

        return PathwayResult(
            ks_stat_b0=ks0,   perm_pval_b0=ksp0,
            cohend_b0=cd0,    cohend_pval_b0=cdp0,  mean_diff_b0=md0,
            ks_stat_b1=ks1,   perm_pval_b1=ksp1,
            cohend_b1=cd1,    cohend_pval_b1=cdp1,  mean_diff_b1=md1,
            significant=sig0 and sig1,
            alpha=self.alpha,
        )

    def _check_fitted(self):
        if not self._fitted:
            raise RuntimeError("Call fit() before accessing results.")
