CanDI.candi module
==================

.. automodule:: CanDI.candi.candi
   :members:
   :undoc-members:
   :show-inheritance:

CanDI.data module
=================
The data class is instantiated at import. This class contains paths to all data downloaded with CanDI.
It has internal methods for loading datasets into memory as pandas dataframes.
There are 3 index tables that candi relies on for fetch all data:

- cell_lines
- genes
- locations

These tables are automatically loaded as pandas dataframes upon import of CanDI
It is highly recommended the user familiarize themself with the columns and indexes of these tables.
All candi classes operate through these index tables.

.. automodule:: CanDI.candi.data
   :members:
   :undoc-members:
   :show-inheritance:
