import sys
import os
from pathlib import Path
from prody import parseDCD, parsePDB, writeNMD, EDA, Ensemble

if len(sys.argv) < 3:
    print(f"USAGE: python {os.path.basename(__file__)} <protein_with_traj.pdb> <output.nmd>")
    sys.exit(1)

pdb_file    = sys.argv[1]
output_file = sys.argv[2]

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
eda_ensemble.calcModes(n_modes=3)

writeNMD(output_file, eda_ensemble[:3], structure)
