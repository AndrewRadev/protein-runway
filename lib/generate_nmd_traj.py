import sys
import os
import itertools
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader as MDAMemoryReader
import numpy as np


class NormalModes:
    def __init__(self, nmd_file, output_file, frame_count=100, magnitude_scale=0.1, mode_count=1):
        """
        NormalModes object with input parameters.
        """
        self.nmd_file = nmd_file
        self.output_file = output_file
        self.frame_count = frame_count
        self.magnitude_scale = magnitude_scale
        self.mode_count = mode_count

        self.coordinates = None
        self.atomnames = None
        self.resnames = None
        self.resids = None
        self.modes = []
        self.trajectory = []

    def parse_nmd_file(self):
        """
        Parse the NMD file and extract atomnames, resnames, resids, coordinates, and modes.
        """
        with open(self.nmd_file) as f:
            for line in f:
                line = line.strip()
                section, _, line = line.partition(' ')

                if section == 'atomnames':
                    self.atomnames = np.array(line.split(' '))

                elif section == 'resnames':
                    self.resnames = np.array(line.split(' '))

                elif section == 'resids':
                    self.resids = np.array(line.split(' '))

                elif section == 'coordinates':
                    coords = (float(c) for c in line.split(' '))
                    self.coordinates = np.array(list(self._batched(coords, n=3, strict=True)))

                elif section == 'mode':
                    mode, _, line = line.partition(' ')
                    magnitude, _, line = line.partition(' ')
                    vector_coordinates = line.split(' ')

                    magnitude = float(magnitude) * self.magnitude_scale
                    vector_coordinates = (magnitude * float(vc) for vc in vector_coordinates)
                    vectors = np.array(list(self._batched(vector_coordinates, n=3, strict=True)))

                    self.modes.append((mode, vectors))

    def validate_modes(self):
        """
        Ensure that the modes' shapes match the coordinates.
        """
        for _, vectors in self.modes:
            assert vectors.shape == self.coordinates.shape

    def generate_trajectory(self):
        """
        Generate the trajectory based on the modes.
        """
        coordinates = self.coordinates.copy()

        # Forward trajectory
        for _ in range(0, self.frame_count // 2):
            for _, vectors in self.modes[:self.mode_count]:
                coordinates = coordinates + vectors
            self.trajectory.append(coordinates)

        # Backward trajectory
        for _ in range(0, self.frame_count // 2):
            for _, vectors in self.modes[:self.mode_count]:
                coordinates = coordinates - vectors
            self.trajectory.append(coordinates)

    def write_trajectory(self):
        """
        Write the trajectory to a PDB file.
        """
        n_atoms = self.coordinates.shape[0]
        u = mda.Universe.empty(
            n_atoms=n_atoms,
            n_residues=n_atoms,
            n_frames=self.frame_count,
            atom_resindex=np.arange(n_atoms),
            trajectory=True,
        )

        u.add_TopologyAttr('names', self.atomnames)
        u.add_TopologyAttr('resids', self.resids)
        u.add_TopologyAttr('resnames', self.resnames)
        u.load_new(np.array(self.trajectory), format=MDAMemoryReader)
        u.atoms.write(self.output_file, frames='all')

    @staticmethod
    def _batched(iterable, n, *, strict=False):
        """
        Batch elements from the iterable into tuples of length n.
        """
        if n < 1:
            raise ValueError('n must be at least one')

        iterator = iter(iterable)

        while batch := tuple(itertools.islice(iterator, n)):
            if strict and len(batch) != n:
                raise ValueError('batched(): incomplete batch')
            yield batch


