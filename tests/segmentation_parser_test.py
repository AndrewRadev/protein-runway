import sys
import os
import unittest

from lib.segmentation_parsers import *

class TestSegmentationParser(unittest.TestCase):

    def setUp(self):
        chainsaw_path = 'test_data/1a3w_noPTM.tsv'
        geostas_path = 'test_data/'
        self.chainsaw = ChainsawParser(chainsaw_path)
        self.geostas = GeostasParser(geostas_path)
        self.assertEqual(self.chainsaw.path, chainsaw_path)
        self.assertEqual(self.geostas.path, geostas_path)

    def test_geostas_chopping(self):
        atom_groups = [[1, 2, 3], [10, 11, 20, 21]]
        chops = self.geostas.generate_geostas_chopping(atom_groups)
        self.assertEqual(chops, '1-3,10-11_20-21')
        self.assertEqual(len(atom_groups) - 1, chops.count(','))

    def test_parsers(self):
        test = self.chainsaw.parse()
        exp_result = [('Chainsaw', '8', '2-84_186-490,90-184,501-581_681-984,587-679,1011-1074_1180-1492,1075-1176,1502-1572_1672-1982,1573-1671')]
        self.assertEqual(test, exp_result)
        self.assertTrue(test[0][1].isdigit())
        self.assertEqual(test[0][0], 'Chainsaw')

        test2 = self.geostas.parse()
        self.assertTrue('GeoStaS' in test2[0][0])
        self.assertEqual(len(test2), len(list(Path(self.geostas.path).glob('clustering_*.json'))))


if __name__ == '__main__':
    unittest.main()
