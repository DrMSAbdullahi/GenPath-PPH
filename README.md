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
Custom plotting functions were created using the above libraries to generate visualizations:
Custom visualization scripts were implemented to illustrate topological features:
- **Persistence Diagrams and Barcodes**: These were plotted using the `plot_persistence_diagram` and `plot_persistence_barcode` functions from the 'Gudhi' library in Python.

### 6. **Permutation Test**  
Permutation tests were performed using Python (version 3.12.6), and the **Benjamini–Hochberg** correction method for multiple testing was applied using the **`statsmodels`** package (version 0.13.5). The False Discovery Rate (FDR) threshold was set to $< 0.05$ throughout.

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

Import the package in Python:

```bash
from genpath_pph import core, utils

# Example: generate allowed paths
paths = utils.generate_allowed_paths(network, max_length=3)

# Example: compute boundaries
boundaries = core.compute_boundaries(paths)
```
Run example scripts:

```bash
python examples/example1.py
python examples/example2.py

These scripts demonstrate PPH calculations, Betti number computation, and visualization for selected pathways.

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

If you use this codebook in your research or publication, please cite:

[Abdullahi, M.S.; Piro, R.M.; Suratanee, A.; Plaimas, K. "GenPath-PPH: Integrating Gene Expression and Pathway Networks via Persistent Path Homology Enhances Detection of Disease-Relevant Pathways". Computational and Structural Biotechnology Journal (Submitted), 2025]

For questions or additional permissions, contact [abdullahi.sirajo@udusok.edu.ng].

---

#### References

[1] Han, Z., Feng, W., Hu, R., Ge, Q., Ma, W., Zhang, W., Xu, S., Zhan, B., Zhang, L., Sun, X. et al. (2021). 'RNA-seq profiling reveals PBMC RNA as a potential biomarker for hepatocellular carcinoma.' Sci. Reports, 11, 17797.

[2] Chen, D., Liu, J., Wu, J., Wei, G.-W., Pan, F., Yau, S.-T. (2023). 'Path topology in molecular and materials sciences'. The Journal of Physical Chemistry Letters 14 (4), 954–964.

[3] Chowdhury, S., Mémoli, F. (2018). 'Persistent path homology of directed networks'. In: Proceedings of the Twenty-Ninth Annual ACM-SIAM Symposium on Discrete Algorithms, pp. 1152–1169.

