"""
Code that encapsulates interactions with an  MDAnalysis universe
"""

from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader as MDAMemoryReader
import numpy as np


class Trajectory:
    """
    A Trajectory object is a wrapper for the MDAnalysis ``Universe`` class. It
    encapsulates some MDA-specific function calls to make it simpler to use for
    the purposes of this codebase.

    It's possible to construct a Trajectory from any MDAnalysis Universe.
    However, it would usually happen from topology and trajectory paths (or a
    single combined PDB). It's also possible to create one from an array of
    alpha-carbon coordinates.

    If needed, more delegation methods can be added from a Trajectory to the
    inner Universe to avoid other classes directly interacting with it.
    """

    @staticmethod
    def from_paths(topology_path, trajectory_path=None):
        """
        Build a :class:`Trajectory` object from a topology file and a
        trajectory file. Mirrors the MDAnalysis Universe construction method.
        """

        if trajectory_path:
            mda_universe = mda.Universe(topology_path, trajectory_path)
        else:
            mda_universe = mda.Universe(topology_path)

        return Trajectory(mda_universe)

    @staticmethod
    def from_ca_frames(data, topology_attr={}):
        """
        A shortcut to create an in-memory MDAnalysis Universe with the given
        data as a list of coordinate matrices, each a list of triplets.

        Each coordinate is considered to be an alpha carbon, so the number of
        residues is set to the number of given atoms. Some additional
        (optional) topology attributes you might provide:

        * names: atom names that match the coordinates (default to "CA")
        * resids: residue ids that correspond to the given atoms
        * resnames: residue names that correspond to the given residue ids

        All other topology attributes are passed along to the MDAnalysis
        Universe method ``add_TopologyAttr``.
        """

        frame_count = len(data)
        n_atoms     = len(data[0])

        u = mda.Universe.empty(
            n_atoms=n_atoms,
            n_residues=n_atoms,
            n_frames=frame_count,
            atom_resindex=np.arange(n_atoms),
            trajectory=True,
        )

        if 'names' not in topology_attr:
            topology_attr['names'] = ['CA'] * n_atoms

        for key, value in topology_attr.items():
            u.add_TopologyAttr(key, value)

        u.load_new(np.array(data), format=MDAMemoryReader)

        return Trajectory(u)

    def __init__(self, mda_universe):
        self.mda_universe = mda_universe

    @property
    def coordinates(self):
        return self.mda_universe.atoms.positions

    @property
    def frames(self):
        return self.mda_universe.trajectory

    def __next__(self):
        """
        Delegates to MDAnalysis to shift the internal trajectory of the
        Universe forward.
        """
        return next(self.mda_universe.trajectory)

    def select_atoms(self, *args):
        """
        Delegates to MDAnalysis to return an AtomGroup.
        """
        return self.mda_universe.select_atoms(*args)

    def write_static(self, path: Path|str, selection='all'):
        """
        Write the static coordinates of the current frame into the given file.
        It's expected that it's a PDB, but anything that MDAnalysis accepts
        will work.
        """
        atoms = self.select_atoms(selection)
        atoms.write(str(path))

    def write_frames(self, path: Path|str, selection='all'):
        """
        Write the full trajectory to a file. It's expected that it's a PDB, but
        anything that MDAnalysis accepts will work.
        """
        atoms = self.select_atoms(selection)
        atoms.write(str(path), frames='all')
