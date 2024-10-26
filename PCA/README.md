## Essential Dynamics Analysis using ProDy

### Setup on the VSC

1. Follow instructions on this page to install miniconda (or any conda version you want) and create a conda environment anywhere in your $VSC_DATA folder: [conda env on the VSC](https://docs.vscentrum.be/software/python_package_management.html#install-python-packages-using-conda)

2. Install necessary packages
```bash
module load Python/3.7.0-foss-2018a
conda install conda-forge::prody
conda install conda-forge::mdtraj
```
Also install numpy>1.2 if it's not there.

3. Run the EDA/PCA in an interactive session. You'll need the pdb structure and a trajectory file (any format handled by mdtraj.mdconvert is fine).
```bash
conda activate <your env name>
# you can check what has been installed in your env with "pip list"
python3 pca.py <input trajectory file>
```

Resources:
* Script is based on this example: [Essential Dynamics Analysis with ProDy](http://www.bahargroup.org/prody/tutorials/trajectory_analysis/eda.html#large-files)
* [How to handle NMD files in VMD](http://prody.csb.pitt.edu/tutorials/nmwiz_tutorial/nmwiz.html)