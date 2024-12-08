"""
Generating normal modes, reading and writing mode data from/to files.
"""

import numpy as np
from pathlib import Path
from prody import parsePDB, writeNMD, EDA, Ensemble

from lib.trajectory import Trajectory
import lib.util as util


def generate_nmd_from_pdb(pdb_path: Path|str, nmd_path: Path|str, mode_count=10):
    """
    Parses a PDB file with a structure and a trajectory and calculates the normal
    modes of the trajectory. Writes them to an NMD file.

    This can be an expensive process, which is why it's encapsulated in a
    function that goes from file to file and doesn't produce an in-memory
    structure. To work with the resulting normal modes, use
    `NormalModes.from_nmd`. which is a cheap operation.
    """
    pdb_path = str(pdb_path)
    nmd_path = str(nmd_path)

    protein_name = Path(pdb_path).stem

    # Limit to alpha carbons to keep a low memory profile
    structure = parsePDB(pdb_path).select('calpha')

    ensemble = Ensemble('%s Structure' % protein_name)

    ensemble.addCoordset(structure)
    ensemble.setCoords(structure)
    ensemble.setAtoms(structure)
    ensemble.superpose()

    eda_ensemble = EDA('%s EDA' % protein_name)
    eda_ensemble.buildCovariance(ensemble)
    eda_ensemble.calcModes(n_modes=mode_count)

    # Pass atoms to writeNMD
    writeNMD(nmd_path, eda_ensemble[:mode_count], structure)


class ValidationError(Exception):
    pass


class NormalModes:
    @staticmethod
    def from_nmd(nmd_path: Path|str):
        nm = NormalModes()
        nm.parse_nmd_file(nmd_path)

        return nm

    def __init__(self):
        """
        NormalModes object with input parameters.
        """
        self.coordinates = None
        self.atomnames = None
        self.resnames = None
        self.resids = None
        self.modes = []

    def parse_nmd_file(self, nmd_path):
        """
        Parse the NMD file and extract atomnames, resnames, resids,
        coordinates, and modes. Structure:

            atomnames CA CA CA ...
            resnames SER ARG LEU ...
            resids 0 1 2 3 4 5 6 7 8 11 12 ... <- Note: possible to have gaps
            coordinates 54.260 50.940 73.060 ...
            mode 1 29.61 -0.008 -0.005 0.012 ...
            mode 2 18.28 -0.009 0.004 -0.008 ...
            mode 3 17.50 -0.027 -0.025 0.023 ...
            mode <N> <magnitude> <x1> <y1> <z1> <x2> <y2> <z2> ...

        There may be additional lines that we currently ignore.
        """
        with open(nmd_path) as f:
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
                    self.coordinates = np.array(self._group_in_threes(coords))

                elif section == 'mode':
                    mode,      _, line = line.partition(' ')
                    magnitude, _, line = line.partition(' ')
                    vector_coordinates = line.split(' ')

                    # Note: magnitude recorded in the file is currently unused,
                    # it sometimes creates trajectories that are too large

                    vector_coordinates = (float(vc) for vc in vector_coordinates)
                    vectors = np.array(self._group_in_threes(vector_coordinates))

                    self.modes.append((mode, vectors))

        self._validate_modes()

    def generate_trajectory(self, frame_count=100, vector_scale=2.5, mode_indices=0):
        """
        Generate a trajectory that visualises these modes.

        The frame count is expected to be an even number. The first frame is
        the initial coordinates, then the animation applies the vectors
        forward, then back. The last frame should also be the initial
        coordinates, provided there are no numerical errors.

        By default, the first mode is used, but you can provide the index in
        `self.modes` to use by changing `mode_indices`. If that input is given
        as a range, the whole range of modes will be used to generate the
        trajectory.

        The `vector_scale` parameter simply controls how much to scale the
        vectors to produce a clearly visible trajectory. Setting it too high
        might result in very large changes to the visualization.
        """
        coordinates = self.coordinates.copy()
        trajectory = [coordinates]

        if type(mode_indices) is range:
            target_modes = self.modes[mode_indices.start:mode_indices.stop]
        else:
            target_modes = [self.modes[mode_indices]]

        # Forward trajectory, half the frames without the first one and the
        # last one:
        for _ in range(0, (frame_count - 2) // 2):
            for _, vectors in target_modes:
                coordinates = coordinates + vectors * vector_scale
            trajectory.append(coordinates)

        # Backward trajectory, half the frames, ending in the original
        # coordinates:
        for _ in range(0, frame_count // 2):
            for _, vectors in target_modes:
                coordinates = coordinates - vectors * vector_scale
            trajectory.append(coordinates)

        trajectory = Trajectory.from_ca_frames(trajectory, topology_attr={
            'names': self.atomnames,
            'resids': self.resids,
            'resnames': self.resnames,
        })

        return trajectory

    def _validate_modes(self):
        """
        Ensure that the modes' shapes match the coordinates.
        """
        for mode, vectors in self.modes:
            if vectors.shape != self.coordinates.shape:
                message = 'Vectors and coordinates mismatch in mode {}: {} != {}'.format(
                    mode, vectors.shape, self.coordinates.shape
                )
                raise ValidationError(message)

    def _group_in_threes(self, flat_coordinates):
        return list(util.batched(flat_coordinates, n=3, strict=True))
