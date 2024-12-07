# NMD structure:
#
# atomnames CA CA CA ...
# resnames SER ARG LEU ...
# resids 0 1 2 3 4 5 6 7 8 11 12 ... <- Note: possible to have gaps
# coordinates 54.260 50.940 73.060 ...
# mode 1 29.61 -0.008 -0.005 0.012 ...
# mode 2 18.28 -0.009 0.004 -0.008 ...
# mode 3 17.50 -0.027 -0.025 0.023 ...
# mode <N> <magnitude> <x1> <y1> <z1> <x2> <y2> <z2> ...

import unittest
import tempfile
import warnings
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from lib.normal_modes import NormalModes, ValidationError
from lib.trajectory import TrajectoryWriter

class TestNormalModes(unittest.TestCase):
    def setUp(self):
        self.root_dir = tempfile.TemporaryDirectory()
        self.root_path = Path(self.root_dir.name)

        warnings.simplefilter("ignore")

    def tearDown(self):
        self.root_dir.cleanup()
        warnings.resetwarnings()

    def test_generate_nmd_from_pdb(self):
        nm = NormalModes()

        # Generate a PDB file from our example inputs for testing:
        pdb_path = self.root_path / 'input.pdb'
        test_top = '01_input/top/5vde_example_complex.top'
        test_traj = '01_input/traj/5vde_example_10-20ns_100snap.trr'
        mda_universe = mda.Universe(test_top, test_traj)
        TrajectoryWriter(mda_universe).write_trajectory_file(pdb_path)

        # Generate normal modes from this input
        nmd_file = self.root_path / 'output.nmd'
        nm.generate_nmd_from_pdb(pdb_path, nmd_file)

        expected_calpha = mda_universe.select_atoms('name = CA')
        self.assertEqual(len(nm.structure.getCoords()), len(expected_calpha))

    def test_parse_nmd_file(self):
        nm = NormalModes()

        nmd_file = self.root_path / "input.nmd"
        self._create_nmd_file(str(nmd_file), """
            coordinates 1 1 1 -1 -1 -1
            mode 1 20.0 1 2 3 4 5 6
            mode 2 10.0 4 5 6 1 1 1
        """)

        nm.parse_nmd_file(nmd_file)

        self.assertTrue(nm.coordinates.tolist(), np.array([(1, 1, 1), (-1, -1, -1)]).tolist())
        self.assertEqual(nm.modes[0][1].tolist(), [[1, 2, 3], [4, 5, 6]])
        self.assertEqual(nm.modes[1][1].tolist(), [[4, 5, 6], [1, 1, 1]])

    def test_generate_trajectory(self):
        nm = NormalModes()

        nm.coordinates = np.random.rand(10, 3)  # Mock coordinates
        nm.modes = [('mode1', np.random.rand(10, 3))]  # Mock modes

        trajectory = nm.generate_trajectory()
        self.assertEqual(trajectory.trajectory.n_frames, nm.frame_count)

    def test_validate_modes(self):
        nm = NormalModes()

        # Test the _validate_modes method
        nm.coordinates = np.random.rand(10, 3)  # Mock coordinates
        nm.modes = [('mode1', np.random.rand(10, 3))]  # Mock modes

        # Should not raise a validation error
        nm._validate_modes()

        with self.assertRaises(ValidationError):
            # One more coordinate than mode vector:
            nm.coordinates = np.random.rand(11, 3)  # Mock coordinates
            nm.modes = [('mode1', np.random.rand(10, 3))]  # Mock modes

            nm._validate_modes()

    def _create_nmd_file(self, name, data):
        with open(self.root_path / name, 'w') as f:
            for line in data.split("\n"):
                f.write(line.strip())
                f.write("\n")


if __name__ == '__main__':
    unittest.main()
