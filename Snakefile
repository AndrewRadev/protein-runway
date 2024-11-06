from snakemake.io import glob_wildcards

pdb_files = glob_wildcards("data/pdb/{pdb_file}.pdb").pdb_file

MERIZO_OK        = 'vendor/merizo_ok.txt'
CHAINSAW_OK      = 'vendor/chainsaw_ok.txt'
MDTASK_OK_PREFIX = 'vendor/mdtask_ok'
MDTASK_OK        = f'{MDTASK_OK_PREFIX}.txt'

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

rule mdtask_setup:
    output: directory('vendor/mdtask')
    shell:
        """
        git clone https://github.com/RUBi-ZA/MD-TASK vendor/mdtask
        cd vendor/mdtask

        # Check out known commit:
        git checkout 6cff460

        # Delete git history to save space:
        rm -rf .git/
        """

rule mdtask_qa:
    input: 'vendor/mdtask'
    output: MDTASK_OK
    shell:
        """
        python vendor/mdtask/calc_correlation.py \
            --trajectory vendor/mdtask/example/example_small.dcd \
            --topology vendor/mdtask/example/example_small.pdb \
            --prefix {MDTASK_OK_PREFIX} \
            --step 100 \
            --lazy-load
        """

rule get_chainsaw_clustering:
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

rule get_merizo_clustering:
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

rule get_mdtask_clustering:
    input:
        correlation_file="output/correlation/{pdb_file}.txt"
    output:
        output_file="output/correlation/{pdb_file}.tsv"
    shell:
        """
        python scripts/segment_by_correlation.py {input.correlation_file} {output.output_file}
        """

rule rename_top_file_to_prmtop:
    input: "{path}.top"
    output: "{path}.prmtop"
    shell: "mv {input} {output}"

rule build_pdb_with_traj:
    input:
        topology="data/top/{protein_name}_complex.prmtop",
        trajectory="data/traj/{protein_name}_10-20ns_100snap.trr"
    output:
        pdb_file="data/pdb/{protein_name}.with_traj.pdb"
    shell: """
        mdconvert {input.trajectory} -t {input.topology} -o {output.pdb_file}
    """

rule build_nmd_file:
    input:
        pdb_file="data/pdb/{protein_name}.with_traj.pdb",
    output:
        nmd_file="output/pca/{protein_name}.nmd"
    shell:
        """
        python scripts/build_nmd_file.py {input.pdb_file} {output.nmd_file}
        """

rule get_trajectory_correlation:
    input:
        topology="data/top/{protein_name}_complex.prmtop",
        trajectory="data/traj/{protein_name}_10-20ns_100snap.trr"
    output:
        correlation_png="output/correlation/{protein_name}.png",
        correlation_txt="output/correlation/{protein_name}.txt"
    shell:
        """
        python vendor/mdtask/calc_correlation.py \
            --trajectory {input.trajectory} \
            --topology {input.topology} \
            --prefix output/correlation/{wildcards.protein_name} \
            --step 100 \
            --lazy-load
        """
