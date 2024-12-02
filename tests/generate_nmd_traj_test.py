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
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from lib.normal_modes import NormalModes
from lib.trajectory import TrajectoryWriter

class TestNormalModes(unittest.TestCase):

    def setUp(self):
        self.nm = NormalModes()

        self.root_dir = tempfile.TemporaryDirectory()
        self.root_path = Path(self.root_dir.name)

    def tearDown(self):
        self.root_dir.cleanup()

    def test_generate_nmd_from_pdb(self):
        # Generate a PDB file from our example inputs for testing:
        pdb_path = self.root_path / 'input.pdb'
        test_top = '01_input/top/5vde_example_complex.top'
        test_traj = '01_input/traj/5vde_example_10-20ns_100snap.trr'
        mda_universe = mda.Universe(test_top, test_traj)
        TrajectoryWriter(mda_universe).write_trajectory_file(pdb_path)

        # Generate normal modes from this input
        nmd_file = self.root_path / 'output.nmd'
        self.nm.generate_nmd_from_pdb(pdb_path, nmd_file)

        expected_calpha = mda_universe.select_atoms('name = CA')
        self.assertEqual(len(self.nm.structure.getCoords()), len(expected_calpha))

    def test_parse_nmd_file(self):
        nm = NormalModes(magnitude_scale=1)

        nmd_file = self.root_path / "input.nmd"
        self._create_nmd_file(str(nmd_file), """
            coordinates 1 1 1 -1 -1 -1
            mode 1 20.0 1 2 3 4 5 6
            mode 2 10.0 4 5 6 1 1 1
        """)

        nm.parse_nmd_file(nmd_file)

        self.assertTrue(nm.coordinates.tolist(), np.array([(1, 1, 1), (-1, -1, -1)]).tolist())
        self.assertEqual(nm.modes[0][1].tolist(), [[20, 40, 60], [80, 100, 120]])
        self.assertEqual(nm.modes[1][1].tolist(), [[40, 50, 60], [10, 10, 10]])

        nm = NormalModes(magnitude_scale=0.1)
        nm.parse_nmd_file(nmd_file)

        self.assertTrue(nm.coordinates.tolist(), np.array([(1, 1, 1), (-1, -1, -1)]).tolist())
        self.assertEqual(nm.modes[0][1].tolist(), [[2, 4, 6], [8, 10, 12]])
        self.assertEqual(nm.modes[1][1].tolist(), [[4, 5, 6], [1, 1, 1]])

    def test_generate_trajectory(self):
        self.nm.coordinates = np.random.rand(10, 3)  # Mock coordinates
        self.nm.modes = [('mode1', np.random.rand(10, 3))]  # Mock modes

        trajectory = self.nm.generate_trajectory()
        self.assertEqual(trajectory.trajectory.n_frames, self.nm.frame_count)

    def test_validate_modes(self):
        # Test the _validate_modes method
        self.nm.coordinates = np.random.rand(10, 3)  # Mock coordinates
        self.nm.modes = [('mode1', np.random.rand(10, 3))]  # Mock modes
        try:
            self.nm._validate_modes()  # Should not raise an assertion error
        except AssertionError:
            self.fail("Validation of modes failed unexpectedly!")

    def _create_nmd_file(self, name, data):
        with open(self.root_path / name, 'w') as f:
            for line in data.split("\n"):
                f.write(line.strip())
                f.write("\n")


if __name__ == '__main__':
    unittest.main()
