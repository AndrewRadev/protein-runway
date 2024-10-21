
from prody import *
from pylab import *
import numpy as np


structure = parsePDB("1fuu_noPTM_tleapout.pdb").select('index 1 to 9510')

print(structure.numAtoms()) #9510

ensemble = parseDCD("1fuu_traj.dcd")

ensemble.setCoords(structure)

ensemble.setAtoms(structure)

ensemble.superpose()

eda_ensemble = EDA('1FUU Ensemble')

eda_ensemble.buildCovariance(ensemble)

eda_ensemble.calcModes(n_modes=10)

for mode in eda_ensemble[:4]:
   print(calcFractVariance(mode).round(2))

writeNMD('mdm2_eda.nmd', eda_ensemble[:3], structure.select('calpha'))

