Usage
=====

To use Protein Runway for visualizing molecular movements, follow these steps:

1. **Prepare Input Data**:
   - Ensure you have trajectory (`.trr`) and topology (`.top`) files.

2. **Run the Workflow**:
   .. code-block:: bash

      snakemake --snakefile

3. **Visualize Results**:
   - Import the 03_output/ `.pdb` files to Blender via Add-on
   - Choose segmentation method (chainsaw, merizo, geostas) and number of clusters

