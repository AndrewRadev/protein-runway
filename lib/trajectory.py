from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader as MDAMemoryReader
import numpy as np


class Trajectory:
    @staticmethod
    def from_paths(topology_path, trajectory_path=None):
        if trajectory_path:
            mda_universe = mda.Universe(topology_path, trajectory_path)
        else:
            mda_universe = mda.Universe(topology_path)

        return Trajectory(mda_universe)

    @staticmethod
    def from_ca_frames(data, topology_attr={}):
        frame_count = len(data)
        n_atoms     = len(data[0])

        u = mda.Universe.empty(
            n_atoms=n_atoms,
            n_residues=n_atoms,
            n_frames=frame_count,
            atom_resindex=np.arange(n_atoms),
            trajectory=True,
        )

        for key, value in topology_attr.items():
            u.add_TopologyAttr(key, value)

        u.load_new(np.array(data), format=MDAMemoryReader)

        return Trajectory(u)

    def __init__(self, mda_universe):
        self.mda_universe = mda_universe

    @property
    def frames(self):
        return self.mda_universe.trajectory

    def select_atoms(self, *args):
        return self.mda_universe.select_atoms(*args)

    def write_static(self, path: Path|str, selection='all'):
        atoms = self.select_atoms(selection)
        atoms.write(str(path))

    def write_frames(self, path: Path|str, selection='all'):
        atoms = self.select_atoms(selection)
        atoms.write(str(path), frames='all')
