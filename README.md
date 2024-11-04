## Setup

This guide will describe how to set up an environment using micromamba, but conda or mamba are also fine. Micromamba seems to be the fastest.

To install micromamba ([Reference](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)):

```bash
curl -L micro.mamba.pm/install.sh > install_micromamba.sh
bash install_micromamba.sh
```

**Important**: When asked for a "root prefix", choose `.micromamba`, without the `~` at the beginning. By default, micromamba creates environments in the home directory. You don't want that on the VSC or you'll run out of quota. Choosing `.micromamba` is going to install packages inside of the repo itself in a folder `.micromamba`, which is gitignored.

Set up and activate a `protein-runway` environment from the env file in the repo:

```bash
micromamba create -f micromamba_env.yml
micromamba activate protein-runway
```

Arguments:

- `-n protein-runway`: Name of the environment
- `-c conda-forge`: Source for installing dependencies
- `python=3.11`: Install a standalone python so we can be independent from the VSC python

Install dependencies:

```bash
pip install -r requirements.txt
```

## Run

Before running any code, we need to activate the micromamba environment:

```bash
micromamba activate protein-runway
```

For now, just a simple test script taken from <https://github.com/MDAnalysis/mdanalysis>

```bash
# Move "crambin_md.pdb" and "crambin_md.xtc" into "data/"
python md_test.py
```
