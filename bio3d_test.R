# Install bio3d via micromamba install r-bio3d
#
# load script in R with:
#
# > source('bio3d_test.R')

library(bio3d)

traj = read.dcd('data/pdb/1a3w_noPTM.with_traj.ca.dcd')
clustering3 = geostas(traj, K=3)

# Reuse computed atomic movement similarity matrix (AMSM):
amsm = clustering3$amsm

clustering4 = geostas(traj, K=4)
clustering5 = geostas(traj, K=5)
clustering6 = geostas(traj, K=6)
