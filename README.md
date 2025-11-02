# GenPath-PPH: Gene ExpressioN and PATHway network integration using Persistent Path Homology

This repository contains the code and supplementary material for the manuscript:
Abdullahi et al., "GenPath-PPH: Integrating Gene Expression and Pathway Networks via Persistent Path Homology Enhances Detection of Disease-Relevant Pathways", submitted to *Computational and Structural Biotechnology Journal*, 2025.

---

## Table of Contents

- [Introduction](#introduction)
- [Datasets](#datasets)
- [GenPath-PPH Analysis](#genpath-pph-analysis)
- [Project Dependencies](#project-dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [License](#license)
- [References](#references)

---

## Introduction

**GenPath-PPH** introduces a novel approach for detecting disease-relevant pathways by integrating gene expression data with biological pathway networks using **persistent path homology (PPH)**. This enables detection of pathways significantly altered in disease conditions, exemplified with **hepatocellular carcinoma (HCC)** from peripheral blood mononuclear cell (PBMC) samples.

---

## Datasets

We used the following datasets in our project. 

- **RNA-Seq Gene Expression**: This comprises 26,575 genes from 34 samples (17 healthy, 17 disease). Original data from Han et al.[1], NCBI project PRJNA739257 [here](https://dataview.ncbi.nlm.nih.gov/object/PRJNA739257). 

- **KEGG Pathways**: We obtained 251 metabolic and signaling pathways from the KEGG database using the packages; 'AnnotationDbi' (v1.66.0) and 'org.Hs.eg.db' (v3.19.1) in R (v4.4.1).

- **Pathways Networks**: Interaction networks (e.g., activations, inhibitions) parsed from KEGG Markup Language (KGML) files downloaded from the KEGG website (https://www.genome.jp/kegg/pathway.html).

---

## GenPath-PPH Analysis

The core methodology of **GenPath-PPH** consists of the following steps for identifying pathways that are most likely implicated in disease conditions.

### 1. **Data Preparation**
- **Data Peprocessing**: Transcript Per Million (TPM) normalization is applied to the raw gene expression data, followed by log₂ transformation to account for transcript length, sequencing depth, and to stabilize variance.

- **Pathway-Specific Gene Expression**: Gene expression subsets corresponding to each biological pathway are extracted by mapping genes to their respective pathways. This yields pathway-specific datasets representing gene activity patterns within each pathway.

- **Condition Grouping**: Samples are grouped into experimental conditions (e.g., **HCC** vs. **Control**) to examine how gene expression differences influence pathway topology.

### 2. **Persistent Path Homology (PPH) Computation**
For each pathway and condition, PPH is computed to capture (track) topological changes (variations) in the gene expression network:
   - **Filtration Scale**: Ranges from 0 to 1, where 0 represents the strongest correlation between genes and 1 represents the weakest.
   - **Step Size**: A step of 0.01 yields 101 Betti numbers per pathway and condition across the filtration scale [0,1].

The output comprises two sets of Betti numbers (for **HCC** and **Control**), describing how the network’s topological features evolve with varying correlation strength. This enables a comparative analysis of pathway network alterations between the two conditions.

### 3. **Statistical Significance**

#### a. **Global Differences**:
Persistence landscapes (PLs) were obtained for each pathway output from the PPH computation. The average landscape was then computed for each group (**HCC** and **Control**). The difference between these average landscapes was quantified using the supremum norm ($\|\cdot\|$∞) and the $L_{1,2}$-norm, and the statistical significance of these differences was evaluated through permutation testing.

#### b. **Pathway-Level Differences**:
Kolmogorov–Smirnov (KS) tests and Cohen’s $d$ were applied to quantify the statistical significance of the differences observed between the two conditions (**HCC** and **Control**) for each pathway through permutation testing.

The obtained $p$-values were corrected for multiple testing using the Benjamini–Hochberg (BH) method. The adjusted $p$-values or False Discovery Rate (FDR) threshold was set to $< 0.05$ throughout.

### 4. **PPH Implementation Details**

- **Dissimilarity Matrix**: This was constructed using $1 - |\rho|$, where $\rho$ is the Pearson correlation coefficient between genes.

- **Path Complex Construction**: For each pathway and condition, path complexes are built using gene point clouds obtained from the dissimilarity matrix and the underlying pathway network structure.

### 5. **Visualization**  

Custom and library-based visualization scripts were implemented to illustrate both network-level and topological features.

- **Custom Visualizations**:  
  - **Correlation Heatmaps** showing inter-gene relationships within each pathway.  
  - **Betti Number Evolution Plots** across filtration values.  
  - **PCA Plots** for sample clustering.  
  - **UpSet Plots** for pathway overlap and intersection analysis.  

- **Library-Based Visualizations (Gudhi)**:  
  - **Persistence Diagrams** using `plot_persistence_diagram`.  
  - **Persistence Barcodes** using `plot_persistence_barcode`.  

These visualizations provide both qualitative and quantitative insights into the topological evolution and expression variability of pathways under different biological conditions.

---

## Project Dependencies

The following dependencies are required for running the analysis:

### 1. **Python Libraries**  
The following Python libraries are required:

- **`persim==0.2.1`**: Used for persistence landscape (PL) related computations.
- **`scipy==1.14.1`**: Used for Kolmogorov-Smirnov (KS) test statistics and effect size computations.
- **`statsmodels==0.13.5`**: Used for applying the Benjamini–Hochberg correction method for multiple testing.
- **`biopython==1.76`**: Used for bioinformatics functionalities.
- **`seaborn==0.11.2`**: Used for advanced statistical visualizations.
- **`matplotlib==3.4.3`**: Used for general plotting.
- **`mpl_toolkits==3.1.1`**: Used for 3D plotting.

### 2. **R Libraries**  
The following R packages are required:

- **`annotationDbi==1.66.0`**: Used for database management and gene annotation.
- **`org.Hs.eg.db==3.19.1`**: Homo sapiens genome annotation package.
- **`edgeR==3.42.4`**: Used for data normalization and transformation.

### 3. **GenPath-PPH Core Framework**  
This framework builds on a modified version of [PathHom](https://github.com/WeilabMSU/PathHom), where correlation-based dissimilarity and filtration-level edge interaction tracking functions were added to enhance the PPH computation.

### 4. **External Tools Used for Comparison**
The following packages were used only for benchmarking or comparison:
- **`DESeq2 (R package)`**: Used for differential gene expression (DEG) analysis.
- **`gseapy==v1.1.3`**: Used for Gene Set Enrichment Analysis (GSEA).

---

## Core Module and Enhancements

The main computation is implemented in `core.py` through the GenPathHomology class (an extended version of PathHomology). This module performs Persistent Path Homology (PPH) analysis on pathway networks by integrating:
- **Point clouds** derived from gene expression data.
- **Directed graphs (digraphs)** obtained from KEGG pathway networks.

Integration is achieved by masking the point cloud with the corresponding digraph, combining both gene-level and network-level information into a unified representation.

### Enhancements in GenPath-PPH

Compared to the original PathHom implementation, this version introduces several improvements:

- **Flexible distance metrics**: Users can choose among multiple distance options for computing pairwise relationships between genes:
   - `euclidean` (used in PathHom)
   - `1 - correlation`
   - `1 - |correlation|` (default for GenPath-PPH)

- **Filtration-level edge tracking**: In addition to Betti number distributions (series), the function also returns the set of edges active (present) at each filtration level, allowing detailed inspection of the evolving topological structure.

- **Integration-ready output**: 
For each pathway and condition, the framework produces:
   - Betti number series describing topological changes across filtrations.
   - Filtration-level edge sets tracking the evolution of gene interactions.

These outputs can be directly integrated into downstream statistical, biological, or visualization analyses (e.g., group comparisons using persistence landscapes, pathway-level significance testing, tracking gene-gene connectivity dynamics across filtrations to identify stable or condition-specific interactions, quantifying edge persistence or node centrality changes to reveal rewired subnetworks, visualizing edge evolution heatmaps, persistence barcodes, or module reorganization across filtrations).

---

## Installation

To use this code, you will need Python and several dependencies. You can install them using `pip`:

```bash
# Clone the repository
git clone https://github.com/DrMSAbdullahi/GenPath-PPH.git
cd GenPath-PPH

# Install the necessary packages
pip install -r requirements.txt
```

---

## Usage

### Import the package in Python

```python
from genpath_pph import core, utils

# Example: generate allowed paths
paths = utils.generate_allowed_paths(network, max_length=2)

# Example: compute boundaries
boundaries = core.compute_boundaries(paths)
```

### Run some examples

**1. Toy example test** This is a small simulated gene expression dataset and toy graph used in our manuscript to demonstrate PPH computation:
```bash
python core.py --test
```
Running above will:
- Generate a small simulated gene expression dataset.
- Compute persistent path homology on a toy graph.
- Print Betti numbers for each dimension across the filtration.

**2. Example for real pathways** This demononstrate PPH computation for real biological pathways (p53 signaling pathway and ferroptosis):

```bash
python examples/p53_signaling.py
python examples/ferroptosis.py
```

**3. Jupyter Notebook examples** – A detailed analysis comparing PH vs PPH on a toy example, including visualization of barcodes. And a code for global dofference analysis and pathway identification via pathway level analysis:

```bash
jupyter notebook notebooks/toy_example_ph_vs_pph.ipynb
jupyter notebook notebooks/genpath_pph_global.ipynb
jupyter notebook notebooks/genpath_pph_pathway_level.ipynb
```

These examples demonstrate how to compute PPH (Betti number computation), visualize topological features, and interpret identified pathways.

---

## Project Structure

```text
GenPath-PPH/
│
├── README.md
├── LICENSE
├── requirements.txt        # Python dependencies
├── genpath_pph/            # Core package
│   ├── __init__.py
│   ├── core.py             # Core PPH computations
│   ├── utils.py            # Utility functions
│   └── analysis.py         # Optional analysis scripts
├── examples/               # Example scripts
│   ├── example1.py
│   └── example2.py
└── tests/                  # Optional unit tests
```

---

## License

### Code License

The code in this repository is licensed under the **MIT License**. You are free to use, modify, and distribute the code, provided that the following conditions are met:

- The original copyright notice and the MIT License text must be included in all copies or substantial portions of the code.
- The code is provided "as is," without warranty of any kind.

For the full MIT License text, please see the [LICENSE](./LICENSE.txt) file.

### Data License

The data in this repository is licensed under the **Creative Commons Attribution 4.0 International License (CC BY 4.0)**. You are free to:

- **Share**: Copy and redistribute the material in any medium or format.
- **Adapt**: Remix, transform, and build upon the material for any purpose, including commercial purposes.

Under the following conditions:

- **Attribution**: You must provide appropriate credit, give a link to the license, and indicate if changes were made. You may do this in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
- **No additional restrictions**: You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

The data is provided "as is," without any warranty of fitness for a particular purpose. Users are solely responsible for their use of the data and any consequences thereof.

For the full Creative Commons Attribution 4.0 International License text, please visit the [CC BY 4.0 License](https://creativecommons.org/licenses/by/4.0/).

---

#### References

[1] Han, Z., Feng, W., Hu, R., Ge, Q., Ma, W., Zhang, W., Xu, S., Zhan, B., Zhang, L., Sun, X. et al. (2021). 'RNA-seq profiling reveals PBMC RNA as a potential biomarker for hepatocellular carcinoma.' Sci. Reports, 11, 17797.

[2] Chen, D., Liu, J., Wu, J., Wei, G.-W., Pan, F., Yau, S.-T. (2023). 'Path topology in molecular and materials sciences'. The Journal of Physical Chemistry Letters 14 (4), 954–964.

[3] Chowdhury, S., Mémoli, F. (2018). 'Persistent path homology of directed networks'. In: Proceedings of the Twenty-Ninth Annual ACM-SIAM Symposium on Discrete Algorithms, pp. 1152–1169.

---

#### Citation

If you use this codebook in your research or publication, please cite:

[Abdullahi, M.S.; Piro, R.M.; Suratanee, A.; Plaimas, K. "GenPath-PPH: Integrating Gene Expression and Pathway Networks via Persistent Path Homology Enhances Detection of Disease-Relevant Pathways". Computational and Structural Biotechnology Journal (Submitted), 2025]

For questions or additional permissions, contact [abdullahi.sirajo@udusok.edu.ng].
