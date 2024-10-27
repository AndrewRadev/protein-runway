import sys
import os
from prody import parseDCD, parsePDB, writeNMD

if len(sys.argv) < 4:
    print(f"USAGE: python {os.path.basename(__file__)} <protein.pdb> <protein.dcd> <output.nmd>")
    sys.exit(1)

pdb_file    = sys.argv[1]
dcd_file    = sys.argv[2]
output_file = sys.argv[3]

ensemble = parseDCD(dcd_file)
num_atoms = ensemble.numAtoms()
structure = parsePDB(pdb_file).select('index 1 to %d' % num_atoms)

ensemble.setCoords(structure)
ensemble.setAtoms(structure)
ensemble.superpose()

eda_ensemble = EDA('%s Ensemble' % wildcards.protein_name)
eda_ensemble.buildCovariance(ensemble)
eda_ensemble.calcModes(n_modes=10) # can do more but will take longer
