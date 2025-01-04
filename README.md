# CanDI - A global cancer data integrator

[![PyPI](https://img.shields.io/pypi/v/PyCanDI)](https://pypi.org/project/PyCanDI/)
[![Downloads](https://static.pepy.tech/badge/pycandi)](https://pepy.tech/project/pycandi)
[![Documentation Status](https://readthedocs.org/projects/candi/badge/?version=latest)](https://candi.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.1186/s13073-021-00987-8.svg)](https://doi.org/10.1186/s13073-021-00987-8)
[![Dataverse](https://img.shields.io/badge/Dataverse-10.7910/DVN/JIAT0H-red)](https://doi.org/10.7910/DVN/JIAT0H)

## Installation

CanDI is now available on [PyPI](https://pypi.org/project/PyCanDI/) and
can be installed with pip. Then, a command from CanDI will automatically
download stable datasets from
[Dataverse](https://doi.org/10.7910/DVN/JIAT0H).

``` bash
# Package Installation & Prepare Datasets
pip install PyCanDI && candi-install
```

Downloaded and formatted datasets would organize this way:

``` 
.
├── config.ini # modified after Installation 
├── depmap
│   ├── CCLE_expression.csv
│   ├── CCLE_fusions.csv
│   ├── CCLE_gene_cn.csv
│   ├── CCLE_mutations.csv
│   ├── CCLE_RNAseq_reads.csv
│   ├── CRISPR_gene_dependency.csv
│   ├── CRISPR_gene_effect.csv
│   └── sample_info.csv
├── genes
│   └── gene_info.csv
└── locations
    └── merged_locations.csv
```

**Note**:

:   *Currently, DepMap API is not available for public use. Therefore,
    we are providing the preprocessed datasets for the users based on
    DepMap 21Q4 release. DepMap API will be available in the future to
    download the latest datasets.*

## Usage

### Import CanDI into python

``` python
from CanDI import candi
```

### CanDI Objects

-   `data` : Container for all candi datasets. All access to datasets go
    through data object.
-   `Gene` : Provides cross dataset indexing from the gene perspective.
-   `CellLine` : Provides cross dataset indexing from the cell line
    perspective.
-   `Cancer` : Provides cross dataset indexing by a group of cell lines
    that are all the same tissue.
-   `Organelle`: Provides cross dataset indexing for a group of genes
    whose proteins localize to the same organelle.
-   `CellLineCluster` : Provides cross dataset indexing for a group of
    user defined cell lines.
-   `GeneCluster` : Provides cross dataset indexing for a group of user
    defined genes.

### Demos

| Name | Description |
|------|-------------|
| Getting Started | [Link to notebook](notebooks/get-started.ipynb) |
| *BRCA* Heatmap | [Link to notebook](notebooks/brca_heatmap.ipynb) |
| *KRAS* and *EGFR* Scatter plot | [Link to notebook](notebooks/kras_egfr_scatter.ipynb) |
| CanDI and DESeq2 | [Link to notebook](notebooks/deseq_setup.ipynb) |

## Citation

If you use CanDI in your research, please cite the following paper:

``` bibtex
Yogodzinski C, Arab A, Pritchard JR, Goodarzi H, Gilbert LA. 
A global cancer data integrator reveals principles of synthetic lethality, sex disparity and immunotherapy. 
Genome Med. 2021;13(1):167. Published 2021 Oct 18. doi:10.1186/s13073-021-00987-8
```
