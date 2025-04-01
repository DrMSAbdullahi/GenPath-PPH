# GenPath-PPH: Gene ExpressioN and PATHway network integration using Persistent Path Homology

This is a supplementary material for codes used in the paper:
Abdullahi et al., "GenPath-PPH: Integrating Gene Expression and Pathway Networks via Persistent Path Homology Enhances Detection of Disease-Relevant Pathways", submitted to Journal of Applied and Computational Topology, 2025.

## Table of Contents

- [Introduction](#introduction)
- [Datasets](#datasets)
- [GenPath-PPH Analysis](#GenPath-PPH-Analysis)
- [Permutation Test](#permutation-test)
- [License](#license)
- [Installation](#installation)


## Introduction

This project contains the code and methodology described in the manuscript **GenPath-PPH**, which introduces a novel approach for detecting disease-relevant pathways by integrating gene expression data with biological pathway networks using persistent path homology (PPH). The approach enhances the ability to detect pathways that are significantly altered in disease conditions, such as hepatocellular carcinoma (HCC) from peripheral blood mononuclear cell (PBMC) samples.

## Datasets

We used the following datasets in our project. 

- **RNA-Seq Gene Expression Dataset**: This dataset comprises 26,575 gene expressions from 34 samples, with 17 representing healthy samples and 17 representing disease samples. The original dataset was sourced from the study by Han et al.[1], titled 'RNA-seq profiling reveals PBMC RNA as a potential biomarker for hepatocellular carcinoma.' To access the original sequencing data, please refer to NCBI project PRJNA739257 [here](https://dataview.ncbi.nlm.nih.gov/object/PRJNA739257).
This project 

- **KEGG Pathways Dataset**: We obtained KEGG pathway datasets from the KEGG database using the R software packages AnnotationDbi (v1.66.0) and org.Hs.eg.db (v3.19.1) with (R v4.4.1). Our focus was specifically on the 251 metabolic and signaling pathways.

- **KEGG Pathways Networks**: Pathway interactions/networks (e.g., activations, inhibitions) were parsed from KEGG Markup Language (KGML) files downloaded from the KEGG website (https://www.genome.jp/kegg/pathway.html).

## GenPath-PPH Analysis

The core methodology consists of the following steps to identify pathways that are most likely to be implicated in diseases.

### 1. **Data Preparation**
- **Data Peprocessing**: The first step involves applying TPM to normalize the gene expression data and log2-transform it to account for transcript length and sequencing depth, and to stabilize variance, respectively.

- **Pathway-Specific Gene Expression Datasets**: Next, the gene expression data specific to each biological pathway is extracted. This is done by identifying the genes involved in each pathway, and then creating pathway-specific datasets based on their expression profiles.

- **Condition Grouping**: The dataset is then divided into different condition groups (e.g., **HCC vs Control**) based on sample conditions. This allows the analysis of how gene expression changes between the two conditions and impacts the pathway structure.

### 2. **Persistent Path Homology (PPH) Calculation**
   For each pathway and each condition (i.e., HCC and Control), PPH is computed to track the topological changes (in both dimensions 0 and 1) in the gene expression network using
   
   - **Filtration Scale**: From 0 to 1, where 0 represents the highest correlation between genes and 1 represents the weakest correlation.
   - **Step Size**: A step size of 0.01 results in 101 Betti numbers for each pathway and each condition across the filtration scale [0,1].

The core output of the PPH computation is a series of Betti numbers, which describe the topological features of the pathway at different filtration values (gene correlation strength). The output consists of two sets of Betti numbers per pathway—one for **HCC** and one for **Control**—enabling a comparative analysis of pathway network alterations between the two conditions.

### 3. **Statistical Significance**

#### a. **Global Differences**:
Calculate the average persistence landscapes for each group (disease and control). Then, compute the difference between them and test for significance using the supremum (sup) and 1/2-norms, via permutation testing.

#### b. **Pathway-Level Differences**:
Use Kolmogorov-Smirnov (KS) tests and Cohen’s $d$ to compute statistical significance via permutation testing.

### Project Dependencies

- **PathHom**: We modify the 'PathHom' [2] code to perform PPH computations and obtain the series of Betti numbers. You can find the original 'PathHom' code at [https://github.com/WeilabMSU/PathHom](https://github.com/WeilabMSU/PathHom).

### Persistent Path Homology Computation

- **Dissimilarity Distance Matrix**: The persistent homology computations were based on the inter-gene dissimilarity distance matrix, calculated using  $1 - |\rho|$, where $\rho$ is the Pearson correlation coefficient.

- **Path Complexes Construction**: We adopted the VR complexes construction method to create the simplicial complexes for our sample point clouds.

### Topological Descriptors

- **User-Friendly Functions**: To simplify the calculation of topological descriptors, we have provided user-friendly functions that assist with the computations and further analysis.

### Visualization

- **Persistence Diagrams and Barcodes**: The persistence diagrams and persistence barcodes were plotted using the 'plot_persistence_diagram' and 'plot_persistence_barcode' functions from the 'Gudhi' library in Python.

- Custom plotting functions were created using the `seaborn` and 'matplotlib.pyplot' and 'mpl_toolkits.mplot3d' libraries to generate additional plots.

## Permutation Test

All permutation tests were performed using Python (version 3.12.6) and the Benjamini–Hochberg correction method was applied for multiple testing using 'statsmodels' package (version 0.13.5). The False Discovery Rate (FDR) threshold was set to $< 0.05$ throughout.

## License

## Code License

The code in this repository is licensed under the **MIT License**. You are free to use, modify, and distribute the code, provided that the following conditions are met:

- The original copyright notice and the MIT License text must be included in all copies or substantial portions of the code.
- The code is provided "as is," without warranty of any kind.

For the full MIT License text, please see the [LICENSE](./LICENSE.txt) file.

---

## Data License

The data in this repository is licensed under the **Creative Commons Attribution 4.0 International License (CC BY 4.0)**. You are free to:

- **Share**: Copy and redistribute the material in any medium or format.
- **Adapt**: Remix, transform, and build upon the material for any purpose, including commercial purposes.

Under the following conditions:

- **Attribution**: You must provide appropriate credit, give a link to the license, and indicate if changes were made. You may do this in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
- **No additional restrictions**: You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

The data is provided "as is," without any warranty of fitness for a particular purpose. Users are solely responsible for their use of the data and any consequences thereof.

For the full Creative Commons Attribution 4.0 International License text, please visit the [CC BY 4.0 License](https://creativecommons.org/licenses/by/4.0/).

If you use this codebook in your research or publication, please cite:

[Abdullahi, M.S.; Piro, R.M.; Suratanee, A.; Plaimas, K. "GenPath-PPH: Integrating Gene Expression and Pathway Networks via Persistent Path Homology Enhances Detection of Disease-Relevant Pathways". Journal of Applied and Computational Topology (Submitted), 2025]

For questions or additional permissions, contact [abdullahi.sirajo@udusok.edu.ng].

#### References

[1] Han, Z., Feng, W., Hu, R., Ge, Q., Ma, W., Zhang, W., Xu, S., Zhan, B., Zhang, L., Sun, X. et al. (2021). 'RNA-seq profiling reveals PBMC RNA as a potential biomarker for hepatocellular carcinoma.' Sci. Reports, 11, 17797.

[2] Chen, D., Liu, J., Wu, J., Wei, G.-W., Pan, F., Yau, S.-T. (2023). 'Path topology in molecular and materials sciences'. The Journal of Physical Chemistry Letters 14 (4), 954–964.

[3] Chowdhury, S., Mémoli, F. (2018). 'Persistent path homology of directed networks'. In: Proceedings of the Twenty-Ninth Annual ACM-SIAM Symposium on Discrete Algorithms, pp. 1152–1169.

## Installation

To use this code, you will need Python and several dependencies. You can install them using `pip`:

```bash
# Clone the repository
git clone https://github.com/DrMSAbdullahi/GenPath-PPH.git
cd GenPath-PPH

# Install the necessary packages
pip install -r requirements.txt
