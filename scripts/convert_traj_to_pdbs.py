"""
Input:
- <protein>:        A protein name to use as prefix for the output
- <topology.top>:   A topology file
- <trajectory.trr>: A trajectory
- <output_dir>:     The directory to write outputs to

Output:
- <protein.with_traj.pdb>:    PDB with a trajectory
- <protein.with_traj.ca.pdb>: Same as above, but reduced to only carbon atoms
- <protein.with_traj.ca.dcd>: Same as above, but in a DCD format
- <protein.static.pdb>:       Static snapshot of the first frame
"""

import sys
import os
from pathlib import Path

import MDAnalysis as mda

if len(sys.argv) < 4:
    print(f"USAGE: python {os.path.basename(__file__)} <protein> <topology.top> <trajectory.trr> <output_dir>")
    sys.exit(1)

protein_name    = sys.argv[1]
topology_file   = sys.argv[2]
trajectory_file = sys.argv[3]
output_dir      = sys.argv[4]

output_static_pdb  = Path(output_dir) / f"{protein_name}.static.pdb"
output_traj_pdb    = Path(output_dir) / f"{protein_name}.with_traj.pdb"
output_traj_ca_pdb = Path(output_dir) / f"{protein_name}.with_traj.ca.pdb"
output_traj_ca_dcd = Path(output_dir) / f"{protein_name}.with_traj.ca.dcd"

u = mda.Universe(topology_file, trajectory_file)

atoms = u.select_atoms('protein')
atoms.write(output_static_pdb)
atoms.write(output_traj_pdb, frames='all')

atoms = u.select_atoms('protein and name is CA')
atoms.write(output_traj_ca_pdb, frames='all')
atoms.write(output_traj_ca_dcd, frames='all')
