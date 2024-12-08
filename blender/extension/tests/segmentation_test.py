import unittest

from lib.segmentation import generate_domain_ranges

class TestSegmentation(unittest.TestCase):
    def test_domain_ranges_from_chopping(self):
        domains = generate_domain_ranges('1-2000')
        self.assertEqual(domains, [[range(1, 2000)]])

        domains = generate_domain_ranges('1-100,200-300')
        self.assertEqual(domains, [
            [range(1, 100)],
            [range(200, 300)]
        ])

        domains = generate_domain_ranges('1-100,100-150_250-300,270-280')
        self.assertEqual(domains, [
            [range(1, 100)],
            [range(100, 150), range(250, 300)],
            [range(270, 280)]
        ])

    def test_single_residue_number_in_domain(self):
        domains = generate_domain_ranges('123')
        self.assertEqual(domains, [[range(123, 123)]])

        domains = generate_domain_ranges('1-100,200')
        self.assertEqual(domains, [
            [range(1, 100)],
            [range(200, 200)]
        ])

        domains = generate_domain_ranges('1-100_200_300-400')
        self.assertEqual(domains, [
            [range(1, 100), range(200, 200), range(300, 400)],
        ])


if __name__ == '__main__':
    unittest.main()
