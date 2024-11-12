from snakemake.io import glob_wildcards

protein_names = glob_wildcards("01_input/traj/{protein_name}_10-20ns_100snap.trr").protein_name

MERIZO_OK   = 'vendor/merizo_ok.txt'
CHAINSAW_OK = 'vendor/chainsaw_ok.txt'

rule all:
    input:
        expand("03_output/{protein_name}.segmentation.tsv", protein_name=protein_names),
        expand("03_output/{protein_name}.nmd_traj.pdb", protein_name=protein_names)

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

rule get_chainsaw_clustering:
    input:
        structure_file="02_intermediate/pdb/{protein_name}.static.pdb",
        chainsaw_ok=CHAINSAW_OK
    output:
        output_file="02_intermediate/chainsaw/{protein_name}.tsv"
    shell:
        """
        python vendor/chainsaw/get_predictions.py \
            --structure_file {input.structure_file} \
            --output {output.output_file}

        # Clean up after chainsaw
        if [ -d results ]; then
          rmdir results/
        fi
        """

rule get_merizo_clustering:
    input:
        structure_file="02_intermediate/pdb/{protein_name}.static.pdb",
        merizo_ok=MERIZO_OK
    output:
        output_file="02_intermediate/merizo/{protein_name}.tsv"
    shell:
        """
        python vendor/merizo/predict.py \
            -d cpu \
            -i {input.structure_file} \
            > {output.output_file}
        """

rule build_pdbs:
    input:
        topology="01_input/top/{protein_name}_complex.top",
        trajectory="01_input/traj/{protein_name}_10-20ns_100snap.trr"
    output:
        traj_pdb_file="02_intermediate/pdb/{protein_name}.with_traj.pdb",
        traj_ca_pdb_file="02_intermediate/pdb/{protein_name}.with_traj.ca.pdb",
        traj_ca_dcd_file="02_intermediate/pdb/{protein_name}.with_traj.ca.dcd",
        static_pdb_file="02_intermediate/pdb/{protein_name}.static.pdb"
    shell: """
        python scripts/convert_traj_to_pdbs.py \
            {wildcards.protein_name} {input.topology} {input.trajectory} \
            02_intermediate/pdb/
    """

rule build_nmd_file:
    input:
        pdb_file="02_intermediate/pdb/{protein_name}.with_traj.pdb",
    output:
        nmd_file="02_intermediate/pca/{protein_name}.nmd"
    shell:
        """
        python scripts/build_nmd_file.py {input.pdb_file} {output.nmd_file}
        """

rule build_nmd_trajectory:
    input:
        nmd_file="02_intermediate/pca/{protein_name}.nmd",
    output:
        traj_file="03_output/{protein_name}.nmd_traj.pdb"
    shell:
        """
        python scripts/generate_nmd_traj.py {input.nmd_file} {output.traj_file}
        """

rule segment_by_bio3d_geostas:
    input:
        dcd_file="02_intermediate/pdb/{protein_name}.with_traj.ca.dcd"
    output:
        clustering=directory("02_intermediate/bio3d_geostas/{protein_name}")
    shell:
        """
        Rscript scripts/segment_with_bio3d_geostas.R {input.dcd_file} {output.clustering}
        """

rule collect_segmentation_intermediates:
    input:
        chainsaw="02_intermediate/chainsaw/{protein_name}.tsv",
        merizo="02_intermediate/merizo/{protein_name}.tsv",
        bio3d_geostas="02_intermediate/bio3d_geostas/{protein_name}/"
    output:
        segmentation="03_output/{protein_name}.segmentation.tsv"
    shell:
        """
        python scripts/collect_segmentation_intermediates.py \
            {input.chainsaw} {input.merizo} {input.bio3d_geostas} \
            {output.segmentation}
        """
