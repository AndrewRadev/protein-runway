import MDAnalysis as mda

# Load simulation results with a single line
u = mda.Universe('data/pdb/1fuu_noPTM_complex.top', 'data/traj/1fuu_noPTM_10-20ns_100snap.trr')

print(u)

# Select atoms
ag = u.select_atoms('name OH')
print("Atoms with name OH")
print(ag)

# Atom data made available as Numpy arrays
print("Positions")
print(ag.positions)

# print("Velocities")
# print(ag.velocities)
# Error: MDAnalysis.exceptions.NoDataError: This Timestep has no velocities

# print("Forces")
# print(ag.forces)
# Error: MDAnalysis.exceptions.NoDataError: This Timestep has no forces

# Iterate through trajectories
for ts in u.trajectory:
    print(ag.center_of_mass())
