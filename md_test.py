import MDAnalysis as mda

# Load simulation results with a single line
u = mda.Universe('data/crambin_md.pdb', 'data/crambin_md.xtc')

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
