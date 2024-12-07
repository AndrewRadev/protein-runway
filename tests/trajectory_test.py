import unittest

from lib.trajectory import Trajectory

class TrajectoryTest(unittest.TestCase):
    def test_basic(self):
        coordinates1 = [[1, 2, 3], [4, 5, 6]]
        coordinates2 = [[2, 3, 4], [3, 4, 5]]

        trajectory = Trajectory.from_ca_frames(
            [coordinates1, coordinates2]
        )

        self.assertEqual(len(trajectory.frames), 2)

        self.assertEqual(trajectory.coordinates.tolist(), coordinates1)
        next(trajectory)
        self.assertEqual(trajectory.coordinates.tolist(), coordinates2)


if __name__ == '__main__':
    unittest.main()
