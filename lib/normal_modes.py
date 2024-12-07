import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader as MDAMemoryReader
import numpy as np
from pathlib import Path
from prody import parsePDB, writeNMD, EDA, Ensemble

import lib.util as util

class ValidationError(Exception):
    pass


class NormalModes:
    def __init__(self, frame_count=100, vector_scale=2.5, mode_count=1):
        """
        NormalModes object with input parameters.
        """
        self.frame_count = frame_count
        self.vector_scale = vector_scale
        self.mode_count = mode_count

        self.coordinates = None
        self.atomnames = None
        self.resnames = None
        self.resids = None
        self.modes = []

        self.protein_name = ''
        self.eda_ensemble = None
        self.structure = None

    def generate_nmd_from_pdb(self, pdb_file: Path|str, nmd_file: Path|str):
        """
        Generate the NMD file from the PDB file.
        """
        pdb_file = str(pdb_file)
        nmd_file = str(nmd_file)

        self.protein_name = Path(pdb_file).stem

        # Limit to alpha carbons to keep a low memory profile
        self.structure = parsePDB(pdb_file).select('calpha')

        ensemble = Ensemble('%s Structure' % self.protein_name)

        ensemble.addCoordset(self.structure)
        ensemble.setCoords(self.structure)
        ensemble.setAtoms(self.structure)
        ensemble.superpose()

        self.eda_ensemble = EDA('%s EDA' % self.protein_name)
        self.eda_ensemble.buildCovariance(ensemble)
        self.eda_ensemble.calcModes(n_modes=10)

        # Pass atoms to writeNMD
        writeNMD(nmd_file, self.eda_ensemble[:10], self.structure)

    def parse_nmd_file(self, nmd_file):
        """
        Parse the NMD file and extract atomnames, resnames, resids,
        coordinates, and modes.
        """
        with open(nmd_file) as f:
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
                    mode, _, line = line.partition(' ')
                    magnitude, _, line = line.partition(' ')
                    vector_coordinates = line.split(' ')

                    # Note: magnitude recorded in the file is currently unused,
                    # it sometimes creates trajectories that are too large

                    vector_coordinates = (float(vc) for vc in vector_coordinates)
                    vectors = np.array(self._group_in_threes(vector_coordinates))

                    self.modes.append((mode, vectors))

        self._validate_modes()

    def generate_trajectory(self):
        """
        Generate the trajectory based on the modes.
        """
        coordinates = self.coordinates.copy()
        trajectory = []

        # Forward trajectory
        for _ in range(0, self.frame_count // 2):
            for _, vectors in self.modes[:self.mode_count]:
                coordinates = coordinates + vectors * self.vector_scale
            trajectory.append(coordinates)

        # Backward trajectory
        for _ in range(0, self.frame_count // 2):
            for _, vectors in self.modes[:self.mode_count]:
                coordinates = coordinates - vectors * self.vector_scale
            trajectory.append(coordinates)

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
        u.load_new(np.array(trajectory), format=MDAMemoryReader)

        return u

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
