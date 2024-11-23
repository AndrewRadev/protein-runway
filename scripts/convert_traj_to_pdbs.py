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

from ..lib.trajectory_conversion import TrajectoryConverter

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
converter = TrajectoryConverter(u)

converter.write_static_file(output_static_pdb,   selection='protein')
converter.write_trajectory_file(output_traj_pdb, selection='protein')

converter.write_static_file(output_traj_ca_pdb,     selection='protein and name is CA')
converter.write_trajectory_file(output_traj_ca_dcd, selection='protein and name is CA')
