rule install_merizo:
    output: directory('vendor/merizo')
    shell:
        """
        # Using our own fork of Merizo to fix an issue with histidine residues:
        git clone https://github.com/AndrewRadev/Merizo vendor/merizo
        cd vendor/merizo

        # Check out known commit:
        git checkout d44ae81

        # Delete git history to save space:
        rm -rf .git/
        """

rule check_merizo:
    input: 'vendor/merizo'
    output: 'vendor/merizo_ok.txt'
    shell:
        """
        python vendor/merizo/predict.py -i vendor/merizo/examples/2xdqA.pdb > {output}
        """

rule install_chainsaw:
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

rule check_chainsaw:
    input: 'vendor/chainsaw'
    output: 'vendor/chainsaw_ok.txt'
    shell:
        """
        python vendor/chainsaw/get_predictions.py \
            --structure_file vendor/chainsaw/example_files/AF-A0A1W2PQ64-F1-model_v4.pdb \
            --output {output}
        """
