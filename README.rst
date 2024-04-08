CanDI - A global cancer data integrator
=======================================

|Documentation Status|
|DOI|
|Dataverse|

Package Installation
--------------------

CanDI is now available on `PyPI <https://pypi.org/project/CanDI/>`_ and can be installed with pip:

.. code:: bash
   pip install CanDI

___
For the latest version (development version) install from GitHub:

.. code:: bash
   pip install git+https://github.com/GilbertLabUCSF/CanDI.git


Prepare Datasets
~~~~~~~~~~~~~~~~

The python command from CanDI will automatically download and modify
datasets.

.. code:: bash

   python CanDI/CanDI/setup/install.py

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

Package Usage
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

.. |Documentation Status| image:: https://readthedocs.org/projects/candi/badge/?version=latest
   :target: https://candi.readthedocs.io/en/latest/?badge=latest

.. |DOI| image:: https://zenodo.org/badge/DOI/10.1186/s13073-021-00987-8.svg
   :target: https://doi.org/10.1186/s13073-021-00987-8

.. |Dataverse| image:: https://img.shields.io/badge/Dataverse-10.7910/DVN/JIAT0H-red
  :target: https://doi.org/10.7910/DVN/JIAT0H
