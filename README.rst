CanDI - A global cancer data integrator
=======================================

|PyPI| |Downloads| |Documentation Status| |DOI| |Dataverse|

Installation
------------

CanDI is now available on `PyPI <https://pypi.org/project/PyCanDI/>`_ and can be installed with pip. 
Then, a command from CanDI will automatically download stable datasets from `Dataverse <https://doi.org/10.7910/DVN/JIAT0H>`_.

.. code:: bash

   # Package Installation
   pip install PyCanDI

   # Prepare Datasets
   candi-install

Downloaded and formatted datasets would organize this way:

.. code::

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


**Note**:
   *Currently, DepMap API is not available for public use. Therefore, we are providing the preprocessed datasets for the users
   based on DepMap 21Q4 release. DepMap API will be available in the future to download the latest datasets.*


Usage
-------------

Import CanDI into python
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from CanDI import candi

CanDI Objects
~~~~~~~~~~~~~

-  ``data`` : Container for all candi datasets. All access to datasets
   go through data object.
-  ``Gene`` : Provides cross dataset indexing from the gene perspective.
-  ``CellLine`` : Provides cross dataset indexing from the cell line
   perspective.
-  ``Cancer`` : Provides cross dataset indexing by a group of cell lines
   that are all the same tissue.
-  ``Organelle``: Provides cross dataset indexing for a group of genes
   whose proteins localize to the same organelle.
-  ``CellLineCluster`` : Provides cross dataset indexing for a group of
   user defined cell lines.
-  ``GeneCluster`` : Provides cross dataset indexing for a group of user
   defined genes.

Demos
~~~~~

.. list-table:: Title
   :widths: 25 50
   :header-rows: 1

   * - Name
     - Description
   
   * - [Getting Started](docs/source/get-started.ipynb)
     - ...
   * - [*BRCA* Heatmap](docs/source/brca_heatmap.ipynb)
     - ...
   * - [*KRAS* and *EGFR* Scatter plot](docs/source/kras_egfr_scatter.ipynb)
     - ...
   * - [CanDI and DESeq2](docs/source/deseq_setup.ipynb)
     - ...


Citation
--------

If you use CanDI in your research, please cite the following paper:

.. code:: bibtex

   Yogodzinski C, Arab A, Pritchard JR, Goodarzi H, Gilbert LA. 
   A global cancer data integrator reveals principles of synthetic lethality, sex disparity and immunotherapy. 
   Genome Med. 2021;13(1):167. Published 2021 Oct 18. doi:10.1186/s13073-021-00987-8



.. |PyPI| image:: https://img.shields.io/pypi/v/PyCanDI
   :target: https://pypi.org/project/PyCanDI/
   
.. |Documentation Status| image:: https://readthedocs.org/projects/candi/badge/?version=latest
   :target: https://candi.readthedocs.io/en/latest/?badge=latest

.. |Downloads| image:: https://static.pepy.tech/badge/pycandi
   :target: https://pepy.tech/project/pycandi

.. |DOI| image:: https://zenodo.org/badge/DOI/10.1186/s13073-021-00987-8.svg
   :target: https://doi.org/10.1186/s13073-021-00987-8

.. |Dataverse| image:: https://img.shields.io/badge/Dataverse-10.7910/DVN/JIAT0H-red
  :target: https://doi.org/10.7910/DVN/JIAT0H
