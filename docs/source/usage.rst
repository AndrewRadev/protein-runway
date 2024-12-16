Usage
=====

To use Protein Runway for visualizing molecular movements, follow these steps:

1. **Prepare Input Data**:

   - Ensure you have trajectory (``.trr``) and topology (``.top``) files.
   - The github repository comes with a file in the input directory called ``5vde_example`` that you can use to do a test run.

2. **Run the Workflow**:

   .. code-block:: bash

      snakemake

3. **Visualize Results**:

   - Import the 03_output/ ``*.nmd_traj.pdb`` files to Blender via Add-on
   - Import the 03_output/ ``*.segmentation.tsv`` files to color functional domains
   - Choose segmentation method (Chainsaw, Merizo, GeoStaS) and number of clusters

Visualization Demo
------------------

.. image:: https://img.youtube.com/vi/aGKgV2fPp-o/maxresdefault.jpg
    :alt: IMAGE ALT TEXT HERE
    :target: https://www.youtube.com/watch?v=aGKgV2fPp-o
