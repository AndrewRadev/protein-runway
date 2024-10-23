
from prody import *
from pylab import *
import numpy as np
import os

# activate conda environment
# module load 
#   cluster/genius/interactive
#   Python/3.7.0-foss-2018a
# conda install conda-forge::prody
# conda install conda-forge::mdtraj (for mdconvert)
# make sure numpy is at least v1.2

# parse command line argument
traj_file = sys.argv[1]
print('Trajectory file:', traj_file)
protein_name = traj_file.split('_')[0]
output_file = traj_file[:-4] + '.dcd'
pdb = protein_name+'_noPTM_tleapout.pdb'


# run mdconvert to get dcd format trajectory
os.system('mdconvert %s -o %s' % (traj_file, output_file))

ensemble = parseDCD(output_file)
num_atoms = ensemble.numAtoms()

structure = parsePDB(pdb).select('index 1 to %d' % num_atoms)

ensemble.setCoords(structure)

ensemble.setAtoms(structure)

ensemble.superpose()


eda_ensemble = EDA('%s Ensemble' % protein_name)

eda_ensemble.buildCovariance(ensemble)

eda_ensemble.calcModes(n_modes=10) # can do more but will take longer

for mode in eda_ensemble[:4]:
   print(calcFractVariance(mode).round(2))

nmd_file = protein_name + '_ensemble.nmd'
writeNMD(nmd_file, eda_ensemble[:3], structure)
