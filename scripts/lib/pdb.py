class PDB:
    def __init__(self, static_pdb, traj_pdb, CA_pdb, CA_dcd):
        self.static = static_pdb
        self.traj = traj_pdb
        self.CA = CA_pdb
        self.CA_dcd = CA_dcd
    
    def build_nmd_file(self, nmd_file):
        from trajectory import Trajectory
        from pathlib import Path
        from prody import parseDCD, parsePDB, writeNMD, EDA, Ensemble

        # if len(sys.argv) < 3:
        #     print(f"USAGE: python {os.path.basename(__file__)} <protein_with_traj.pdb> <output.nmd>")
        #     sys.exit(1)

        pdb_file    = self.static
        output_file = nmd_file

        protein_name = Path(pdb_file).stem

        structure = parsePDB(pdb_file)

        # Limit to alpha carbons to keep a low memory profile
        structure = structure.select('calpha')

        ensemble = Ensemble('%s Structure' % protein_name)

        ensemble.addCoordset(structure)
        ensemble.setCoords(structure)
        ensemble.setAtoms(structure)
        ensemble.superpose()

        eda_ensemble = EDA('%s EDA' % protein_name)
        eda_ensemble.buildCovariance(ensemble)
        eda_ensemble.calcModes(n_modes=10)

        writeNMD(output_file, eda_ensemble[:10], structure)

        return Trajectory(nmd_file)