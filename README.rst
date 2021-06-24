CanDI - A global cancer data integrator
=======================================

|Documentation Status|

Package Installation
--------------------

First, you need to clone this repository to use CanDI.

.. code:: bash

   git clone https://github.com/GilbertLabUCSF/CanDI.git

We suggest to use `Conda <https://docs.conda.io/en/latest/>`__ as a
package manager and environment management system. You can create a
fresh conda environment with all ``CanDI``\ ’s requirements using bellow
command:

.. code:: bash

   conda env create -f CaDI/candi.yml -n candi

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

To import ``CanDI``, your active directory in python must be same as the
cloned folder.

.. code:: python

   from CanDI import candi

**OR**, you can add path to the `CanDI` directory if you want to use it from other directories.

.. code:: python

   import sys
   sys.path.append("path-to-candi-directory")

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