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
import numpy as np

from lib.normal_modes import NormalModes

class TestNormalModes(unittest.TestCase):

    def setUp(self):
        self.nm = NormalModes()

    def test_generate_nmd_from_pdb(self):
        pdb_file = 'test_data/1a3w_noPTM.with_traj.pdb'
        nmd_file = 'test_data/1a3w_noPTM.with_traj.nmd'
        self.nm.generate_nmd_from_pdb(pdb_file, nmd_file)
        self.assertIsNotNone(self.nm.structure)
        self.assertIsNotNone(self.nm.eda_ensemble)

    def test_parse_nmd_file(self):
        nmd_file = 'test_data/1a3w_noPTM.with_traj.nmd'
        self.nm.parse_nmd_file(nmd_file)
        self.assertIsNotNone(self.nm.atomnames)
        self.assertIsNotNone(self.nm.resnames)
        self.assertIsNotNone(self.nm.resids)
        self.assertIsNotNone(self.nm.coordinates)
        self.assertGreater(len(self.nm.modes), 0)

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


if __name__ == '__main__':
    unittest.main()
