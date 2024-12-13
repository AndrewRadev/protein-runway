Protein Runway
==============

Protein Runway is a tool for the schematic visualization of large molecules' most important movements based on computational analysis. It consists of:

* A snakemake pipeline that takes a topology and trajectory. It outputs a PDB with a simplified trajectory visualizing its normal modes and a TSV file with various potential segmentations of the protein.
* An object-oriented library that encapsulates the specific processing logic that the workflow applies.
* A blender extension that visualizes the results.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction
   installation
   tutorial
