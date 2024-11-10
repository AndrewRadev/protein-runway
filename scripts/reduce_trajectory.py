"""
Reduce the given trajectory to only its alpha carbons

Input:
- <protein.with_traj.pdb>:    A PDB with an embedded trajectory
- [protein.with_traj.ca.pdb]: The output filename to write, optional, defaults to the input filename with `.ca.`
"""

import sys
import os
from pathlib import Path

import MDAnalysis as mda

if len(sys.argv) < 2:
    print(f"USAGE: python {os.path.basename(__file__)} <protein.with_traj.pdb> [protein.with_traj.ca.pdb]")
    sys.exit(1)

pdb_file    = sys.argv[1]
output_file = sys.argv[2] if len(sys.argv) > 2 else Path(pdb_file).with_suffix('.ca.pdb')

u = mda.Universe(pdb_file)

atoms = u.select_atoms('protein and name is CA')
atoms.write(output_file, frames='all')
