import unittest
import tempfile
import warnings
from pathlib import Path

import numpy as np

from lib.normal_modes import generate_nmd_from_pdb, NormalModes, ValidationError
from lib.trajectory import Trajectory

class NormalModesTest(unittest.TestCase):
    def setUp(self):
        self.root_dir = tempfile.TemporaryDirectory()
        self.root_path = Path(self.root_dir.name)

        warnings.simplefilter("ignore")

    def tearDown(self):
        self.root_dir.cleanup()
        warnings.resetwarnings()

    def test_generate_nmd_from_pdb(self):
        # Generate a PDB file from our example inputs for testing:
        pdb_path = self.root_path / 'input.pdb'
        test_top = '01_input/top/5vde_example_complex.top'
        test_traj = '01_input/traj/5vde_example_10-20ns_100snap.trr'

        trajectory = Trajectory.from_paths(test_top, test_traj)
        trajectory.write_frames(pdb_path)

        # Generate normal modes from this input
        nmd_file = self.root_path / 'output.nmd'
        generate_nmd_from_pdb(pdb_path, nmd_file)

        nm = NormalModes.from_nmd(nmd_file)

        expected_calpha = trajectory.select_atoms('name = CA')
        self.assertEqual(len(nm.coordinates), len(expected_calpha))


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

        trajectory = nm.generate_trajectory(frame_count=10)
        self.assertEqual(len(trajectory.frames), 10)

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
