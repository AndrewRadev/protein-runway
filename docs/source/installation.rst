Setup
=====

This guide will describe how to set up an environment using micromamba, but conda or mamba are also fine. Micromamba seems to be the fastest.

To install micromamba (`Reference <https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html>`_):

.. code-block:: bash
    # If on the VSC:
    cd $VSC_DATA

    curl -L micro.mamba.pm/install.sh > install_micromamba.sh
    bash install_micromamba.sh

This should show you a step-by-step wizard:

.. code-block:: text
    Micromamba binary folder? [~/.local/bin] -> Default is fine
    Init shell (bash)? [Y/n]                 -> Default is fine (may be different based on your default shell)
    Configure conda-forge? [Y/n]             -> Default is fine
    Prefix location? [~/micromamba]          -> Change, see below

**Important**: When asked for a "root prefix", you can choose:

- `.micromamba`, without the `~` at the beginning. This will create one folder named `.micromamba` where you currently are, which contains the default environment. When you create a new environment while inside a project, it will create a new `.micromamba` folder inside that project for your new environment.
- `/some/absolute/path/micromamba`: This will create one directory that will hold all your environments. You can `echo $VSC_DATA/micromamba` and then copy that.

By default, micromamba creates environments in the home directory. You don't want that on the VSC or you'll run out of quota. The `.micromamba` directory is gitignored in the protein-runway project, so it's a fine choice. It's also okay to choose `$VSC_DATA/micromamba` or whatever you like inside of `$VSC_DATA`, just make sure to *expand* it beforehand.

Set up and activate a `protein-runway` environment from the env file in the repo:

.. code-block:: bash
    micromamba create -f micromamba_env.yml
    micromamba activate protein-runway

Arguments:

- `-n protein-runway`: Name of the environment
- `-c conda-forge`: Source for installing dependencies
- `python=3.11`: Install a standalone python so we can be independent from the VSC python

Install dependencies:

.. code-block:: bash
    pip install -r requirements.txt

Annoyingly, this will likely still show some dependency errors, but they don't seem to stop the pipeline from working.

Run snakemake
=============

Before running any code, we need to activate the micromamba environment:

.. code-block:: bash
    micromamba activate protein-runway

The input data consists of two files per protein from Zenodo <https://zenodo.org/records/4650406>_:

- Trajectory file from `Trajectory_snapshots_<code>_<code>.zip`:
    - e.g. `1fuu_noPTM_10-20ns_100snap.trr`
    - Location: `01_input/traj/`
- Topology file from `Topology_stripped.zip`:
    - Example: `1fuu_noPTM_complex.top`
    - Location: `01_input/top/`

At that point, `snakemake -call` should build all necessary output files for these proteins into `03_output`.

Build blender extension
========================

You need to install packages and zip them. To install, this should work from the root of the project:

.. code-block:: bash
    pip wheel mdanalysis -w ./blender/extension/wheels/

Then zipping the extension can be done by running:

.. code-block:: bash
    bash blender/build.sh

But this assumes you're on Linux or Mac, so you can just take the contents of the "blender/extension" folder and zip them manually. You can then install the addon from `Edit > Preferences > Get Extensions`. Then, in the right-hand corner menu, `Install from Disk`.

What blender does then is it just unzips that into a local folder. For me, that folder is `~/.config/blender/4.2/extensions/user_default`. So if you want, you can symlink the extension directly there, so that when you edit the code, you can just re-launch blender. For me, the directory looks like this:

.. code-block:: text
    % ls -al
    total 12K
    drwxr-xr-x 3 andrew andrew 4.0K Nov 12 20:47 .
    drwxr-xr-x 6 andrew andrew 4.0K Nov 12 20:02 ..
    drwxr-xr-x 3 andrew andrew 4.0K Nov 12 20:41 .blender_ext
    lrwxrwxrwx 1 andrew andrew   55 Nov 12 20:47 protein_runway -> /home/andrew/projects/protein-runway/blender/extension/
