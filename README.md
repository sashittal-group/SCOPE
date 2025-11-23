# Exploring the Space of Tumor Phylogenies Consistent with Single-Cell Whole-Genome Sequencing Data
![SCOPE Overview](assets/Overview.png)

**Overview of SCOPE algorithm.** (a) Single-cell whole genome sequencing (scWGS) ultra-low coverage measurement of mutations in individual cells which enables reliable measurement of copy-numbers but complicates measurement of SNV clones. SCOPE takes a two-step approach to infer tumor phylogenies from scWGS data. (b) First, SCOPE estimates cell fractions of mutations within each copy-number cluster using a probabilistic read counts model and identifies forbidden ancestral pairs D of copy-number clusters. (c) Second, SCOPE enumerates all tumor phylogenies consistent with the estimated cell fractions and forbidden ancestral pairs under a copy-number constrained perfect phylogeny model.


## Overview
- Introduces SCOPE, an algorithmic framework that enumerates all tumor phylogenies supported by ultra-low-coverage scWGS data under copy-number–constrained perfect phylogeny assumptions.

- Uses mutation cell-fraction estimates within copy-number clusters and derives necessary and sufficient constraints that these fractions must satisfy to admit a valid phylogeny.

- Demonstrates improved accuracy and runtime over existing tools on simulations and real scWGS datasets, and quantifies phylogenetic uncertainty across cancer samples by identifying when multiple phylogenies are equally supported by the data.

## Installation

```bash
git clone https://github.com/sashittal-group/SCOPE.git
cd SCOPE
pip install .
```

## Run SCOPE on Ovarian Cancer Dataset
This notebook [SCOPE on Ovarian Cancer Dataset](https://github.com/sashittal-group/SCOPE/blob/master/notebooks/Ovarian%20Cancer%20Dataset.ipynb) shows how to run SCOPE on Ovarian Cancer Dataset.
