"""
GenPath-PPH: Gene ExpressioN and PATHway network integration
using Persistent Path Homology.

This package provides tools for:
- Computing Persistent Path Homology (PPH) on directed KEGG pathway graphs.
- Generating Betti number series across expression-induced filtrations.
- Statistical analysis (KS test, Cohen's d, permutation testing).
- Utilities for adjacency extraction and pathway-level summaries.

Reference
---------
Abdullahi, M.S., Piro, R.M., Suratanee, A., Plaimas, K. (2025).
GenPath-PPH: Integrating gene expression and pathway networks via
persistent path homology enhances detection of disease-relevant pathways.
Computational and Structural Biotechnology Journal, 27, 5348-5362.
https://doi.org/10.1016/j.csbj.2025.11.018

Quick start
-----------
    # One-liner: single pathway
    from genpath_pph import run_pathway
    result = run_pathway(X_disease, X_control, adj)
    print(result)

    # One-liner: batch over many pathways
    from genpath_pph import run_batch
    results_df = run_batch(pathway_ids, adj_matrices, pathway_exprs, class_size=17)

    # Object-oriented
    from genpath_pph import GenPathAnalysis
    model = GenPathAnalysis()
    model.fit(X_disease, X_control, adj)
    b0_d, b1_d = model.betti_series("disease")
    result = model.test()

    # Low-level (original classes, unchanged)
    from genpath_pph import GenPathHomology, PathwayDataProcessor, load_toy_data
"""

# ── Core classes (unchanged from original implementation) ─────────────────────
from .core import GenPathHomology, PathwayDataProcessor, load_toy_data

# ── Utility functions (unchanged) ─────────────────────────────────────────────
from .utils import (
    extract_adjacency_matrices,
    extract_pathway_expressions,
    perform_ks_and_effectsize_tests,
    permutation_test_effect_size,
    pooled_std,
    plot_betti_numbers,
    plot_persistent_diagrams,
    plot_persistent_barcode,
    betti_to_pairs_single,
)

# ── Analysis functions (unchanged) ────────────────────────────────────────────
from .analysis import (
    pathway_level_difference_analysis,
    global_difference_analysis,
)

# ── High-level API (new clean interface) ──────────────────────────────────────
from .api import (
    run_pathway,        # One-liner: full pipeline for a single pathway
    run_batch,          # One-liner: full pipeline over many pathways → DataFrame
    GenPathAnalysis,    # Object-oriented interface
    PathwayResult,      # Result dataclass
)

__version__ = "1.0.4"
__author__  = "Muhammad Sirajo Abdullahi"
__email__   = "abdullahi.sirajo@udusok.edu.ng"

__all__ = [
    # Core (original, unchanged)
    "GenPathHomology",
    "PathwayDataProcessor",
    # Utils (original, unchanged)
    "extract_adjacency_matrices",
    "extract_pathway_expressions",
    "perform_ks_and_effectsize_tests",
    "permutation_test_effect_size",
    "pooled_std",
    "plot_betti_numbers",
    "plot_persistent_diagrams",
    "plot_persistent_barcode",
    "betti_to_pairs_single",
    # Analysis (original, unchanged)
    "pathway_level_difference_analysis",
    "global_difference_analysis",
    # High-level API (new)
    "run_pathway",
    "run_batch",
    "GenPathAnalysis",
    "PathwayResult",
]
