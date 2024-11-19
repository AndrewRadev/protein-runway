import unittest

from lib.segmentation import generate_domain_ranges

class TestSegmentation(unittest.TestCase):
    def test_domain_ranges_from_chopping(self):
        domains = generate_domain_ranges('1-2000')
        self.assertEqual(domains, [[range(0, 2000)]])

        domains = generate_domain_ranges('1-100,200-300')
        self.assertEqual(domains, [
            [range(0, 100)],
            [range(199, 300)]
        ])

        domains = generate_domain_ranges('1-100,100-150_250-300,270-280')
        self.assertEqual(domains, [
            [range(0, 100)],
            [range(99, 150), range(249, 300)],
            [range(269, 280)]
        ])


if __name__ == '__main__':
    unittest.main()
