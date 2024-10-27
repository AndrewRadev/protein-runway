from snakemake.io import glob_wildcards

pdb_files = glob_wildcards("data/pdb/{pdb_file}.pdb").pdb_file

MERIZO_OK   = 'vendor/merizo_ok.txt'
CHAINSAW_OK = 'vendor/chainsaw_ok.txt'

rule all:
    input:
        expand("output/merizo/{pdb_file}.tsv", pdb_file=pdb_files),
        expand("output/chainsaw/{pdb_file}.tsv", pdb_file=pdb_files)

rule merizo_setup:
    output: directory('vendor/merizo')
    shell:
        """
        git clone https://github.com/psipred/Merizo vendor/merizo
        cd vendor/merizo

        # Check out known tag:
        git checkout v1.0.0

        # Delete git history to save space:
        rm -rf .git/
        """

rule merizo_qa:
    input: 'vendor/merizo'
    output: MERIZO_OK
    shell:
        """
        python vendor/merizo/predict.py -d cpu -i vendor/merizo/examples/2xdqA.pdb > {output}
        """

rule chainsaw_setup:
    output: directory('vendor/chainsaw')
    shell:
        """
        git clone https://github.com/JudeWells/Chainsaw vendor/chainsaw
        cd vendor/chainsaw

        # Check out known commit:
        git checkout 9ced6e6

        # Delete git history to save space:
        rm -rf .git/

        # Compile stride:
        cd stride
        tar -xzf stride.tgz
        make
        """

rule chainsaw_qa:
    input: 'vendor/chainsaw'
    output: CHAINSAW_OK
    shell:
        """
        python vendor/chainsaw/get_predictions.py \
            --structure_file vendor/chainsaw/example_files/AF-A0A1W2PQ64-F1-model_v4.pdb \
            --output {output}
        """

rule get_chainsaw_predictions:
    input:
        structure_file="data/pdb/{pdb_file}.pdb",
        chainsaw_ok=CHAINSAW_OK
    output:
        output_file="output/chainsaw/{pdb_file}.tsv"
    shell:
        """
        python vendor/chainsaw/get_predictions.py \
            --structure_file {input.structure_file} \
            --output {output.output_file}
        """

rule get_merizo_predictions:
    input:
        structure_file="data/pdb/{pdb_file}.pdb",
        merizo_ok=MERIZO_OK
    output:
        output_file="output/merizo/{pdb_file}.tsv"
    shell:
        """
        python vendor/merizo/predict.py \
            -d cpu \
            -i {input.structure_file} \
            > {output.output_file}
        """

rule build_dcd_file:
    input:
        traj="data/traj/{protein_name}.xtc"
    output:
        dcd_file="output/pca/{protein_name}.dcd"
    shell:
        """
        mdconvert {input.traj} -o {output.dcd_file}
        """

rule build_nmd_file:
    input:
        pdb_file="data/pdb/{protein_name}.pdb",
        dcd_file="output/pca/{protein_name}.dcd"
    output:
        nmd_file="output/pca/{protein_name}.nmd"
    shell:
        """
        python scripts/build_nmd_file.py {input.pdb_file} {input.dcd_file} {output.nmd_file}
        """
