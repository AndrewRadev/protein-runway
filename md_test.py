import MDAnalysis as mda

# Load simulation results with a single line
u = mda.Universe('01_input/top/1fuu_noPTM_complex.top', '01_input/traj/1fuu_noPTM_10-20ns_100snap.trr')
# u = mda.Universe('1a3w_noPTM.nmd_traj.pdb.gz')

print(u)

# Select atoms
ag = u.select_atoms('protein')
print(ag)

# # Atom data made available as Numpy arrays
# print("Positions")
# print(ag.positions)
#
# # Iterate through trajectories
# for ts in u.trajectory:
#     print(ag.center_of_mass())
