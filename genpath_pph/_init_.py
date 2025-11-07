# __init__.py for GenPath-PPH

"""
GenPath-PPH: Gene Expression and PATHway network integration using Persistent Path Homology.

This package provides tools for:
- Computing Persistent Path Homology (PPH) on pathway networks.
- Generating Betti numbers and persistence landscapes.
- Statistical analysis (KS test, Cohen's d, permutation testing) on pathway topologies.
- Utilities for allowed paths, boundaries, and pathway-level summaries.
"""

# Core PPH computations, Betti number calculations
from .core import GenPathHomology, PathwayDataProcessor

# Utilities for path and boundary computations
from .utils import extract_adjacency_matrices, extract_pathway_expressions, compute_betti_numbers

# Optional analysis functions/workflows
from .analysis import pathway_level_difference_analysis, global_difference_analysis

# Define what is exported when using 'from genpath_pph import *'
__all__ = [
    "GenPathHomology",
    "PathwayDataProcessor",
    "extract_adjacency_matrices",
    "extract_pathway_expressions",
    "compute_betti_numbers",
    "utils_generate_allowed_paths",
    "utils_unlimited_boundary_operator",
    "pathway_level_difference_analysis",
    "global_difference_analysis"
]
