[![Tests](https://github.com/AndrewRadev/protein-runway/actions/workflows/tests.yml/badge.svg)](https://github.com/AndrewRadev/protein-runway/actions/workflows/tests.yml)
[![Snakemake workflow](https://github.com/AndrewRadev/protein-runway/actions/workflows/snakemake.yml/badge.svg)](https://github.com/AndrewRadev/protein-runway/actions/workflows/tests.yml)
[![Blender extension](https://github.com/AndrewRadev/protein-runway/actions/workflows/blender_extension.yml/badge.svg)](https://github.com/AndrewRadev/protein-runway/actions/workflows/tests.yml)

To know more about the project:

* 🖼️ See our [poster](./docs/source/resources/poster_2000p.jpeg) overview
* 📄 Read the paper: [PDF](https://raw.githubusercontent.com/AndrewRadev/protein-runway/994d3d357c0a4dbf8e499ebd5a117f12340f154d/docs/source/resources/paper.pdf)
* 📚 Take a look at the documentation: <https://andrewradev.github.io/protein-runway/>
* 📺 Watch a video demonstration: <https://www.youtube.com/watch?v=aGKgV2fPp-o>

The rest of the README consists of installation instructions if you'd like to get it running yourself.

## Setup snakemake workflow

This guide will describe how to set up an environment using micromamba, but conda or mamba are also fine. Micromamba seems to be the fastest. We'll add some specific advice for setting it up on a [VSC account](https://docs.vscentrum.be/index.html).

To install micromamba ([Reference](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)):

```bash
# If on the VSC:
cd $VSC_DATA

# Download installation script
curl -L micro.mamba.pm/install.sh > install_micromamba.sh

# Ideally, skim the script to check it for issues, then:
bash install_micromamba.sh
```

This should show you a step-by-step wizard:

```
Micromamba binary folder? [~/.local/bin]
Init shell (bash)? [Y/n]
Configure conda-forge? [Y/n]
Prefix location? [~/micromamba] -> You might want to change this, see below
```

**Important for the VSC**: When asked for a \"root prefix\", you can choose:

- `.micromamba`, without the `~` at the beginning. This will create one folder named `.micromamba` where you currently are, which contains the default environment. When you create a new environment while inside a project, it will create a new `.micromamba` folder inside that project for your new environment.
- `/some/absolute/path/micromamba`: This will create one directory that will hold all your environments. On the VSC, you can `echo $VSC_DATA/micromamba` and then copy that.

By default, micromamba creates environments in the home directory. You don't want that on the VSC or you'll run out of quota. The `.micromamba` directory is gitignored in the protein-runway project, so it's a fine choice. It's also okay to choose `$VSC_DATA/micromamba` or whatever you like inside of `$VSC_DATA`, just make sure to *expand* it beforehand.

Set up and activate a `protein-runway` environment from the env file in
the repo:

```bash
micromamba create -f micromamba_env.yml
micromamba activate protein-runway
```

Arguments:

-   `-n protein-runway`: Name of the environment
-   `-c conda-forge`: Source for installing dependencies
-   `python=3.11`: The specific version that is required for the project

Install dependencies:

```bash
pip install -r requirements.txt
```

This will likely still show some dependency errors, but they don't seem
to stop the pipeline from working. An alternative setup for the
repository might be to check out the individual projects' requirements
in separate virtual environments or even micromamba environments for
better isolation.

## Run snakemake workflow

Before running any code, we need to activate the micromamba environment:

```bash
micromamba activate protein-runway
```

The input data consists of two files per protein from Zenodo <https://zenodo.org/records/4650406>:

- Trajectory file from `Trajectory_snapshots_<code>_<code>.zip`:
    - e.g. `1fuu_noPTM_10-20ns_100snap.trr`
    - Location: `01_input/traj/`
- Topology file from `Topology_stripped.zip`:
    - Example: `1fuu_noPTM_complex.top`
    - Location: `01_input/top/`

At that point, `snakemake` should build all necessary output files for these proteins into `03_output`.

## Build blender extension

The easiest way to get the blender extension is to use the version that is built by the github actions pipeline, an archive attached to the \"Releases\" of the project: <https://github.com/AndrewRadev/protein-runway/releases>.

If you'd like to set it up yourself from the source, you need to install packages and zip them. From the root of the project:

``` bash
pip wheel MDAnalysis -w ./blender/extension/wheels/
```

Then zipping the extension can be done by running:

```
bash blender/build.sh
```

But this assumes you're on Linux or Mac, so if you're on Windows, you can just take the contents of the "blender/extension" folder and zip them manually.

Either way, you can then install the addon from `Edit > Preferences > Get Extensions`. Then, in the right-hand corner menu, `Install from Disk`.

What blender does is it just unzips that into a local folder. On Linux, that folder is `~/.config/blender/4.2/extensions/user_default`. So if you are actively changing the extension, you can symlink the local folder directly there, so that when you edit the code, you can just re-launch blender.
