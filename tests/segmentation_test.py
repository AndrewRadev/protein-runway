import unittest
import tempfile
import csv
import json
import warnings
from pathlib import Path

import numpy as np

from lib.trajectory import Trajectory
from lib.segmentation import (
    chainsaw,
    merizo,
    geostas,
    ParsingError,
)

class SegmentationTest(unittest.TestCase):
    def setUp(self):
        self.root_dir = tempfile.TemporaryDirectory()
        self.root_path = Path(self.root_dir.name)
        warnings.simplefilter("ignore")

    def tearDown(self):
        self.root_dir.cleanup()
        warnings.resetwarnings()

    def test_chainsaw_parser_success(self):
        self._create_chainsaw_clustering_file(
            'chainsaw.tsv',
            ['_', '_', 6, 2, '1-3,10-12', '_', '_']
        )

        parser = chainsaw.Parser(self.root_path / 'chainsaw.tsv')
        result = list(parser.parse())

        exp_result = [(
            'Chainsaw',
            '2',
            '1-3,10-12',
        )]

        self.assertEqual(result, exp_result)

    def test_chainsaw_parser_failure(self):
        parser = chainsaw.Parser(self.root_path / 'nonexistent.tsv')
        with self.assertRaises(ParsingError):
            list(parser.parse())

        # Mix up the data files:
        self._create_merizo_clustering_file(
            'merizo.tsv',
            ['incorrect data format']
        )

        parser = chainsaw.Parser(self.root_path / 'merizo.tsv')
        with self.assertRaises(ParsingError):
            list(parser.parse())

    def test_merizo_parser_success(self):
        self._create_merizo_clustering_file(
            'merizo.tsv',
            ['_', '_', '_', '_', 2, '_', '_', '1-3,10-12']
        )

        parser = merizo.Parser(self.root_path / 'merizo.tsv')
        result = list(parser.parse())

        exp_result = [('Merizo', '2', '1-3,10-12')]

        self.assertEqual(result, exp_result)

    def test_merizo_parser_failure(self):
        parser = merizo.Parser(self.root_path / 'nonexistent.tsv')
        with self.assertRaises(ParsingError):
            list(parser.parse())

        self._create_merizo_clustering_file(
            'merizo.tsv',
            ['<protein_name>', 'Segmentation failure']
        )

        parser = merizo.Parser(self.root_path / 'merizo.tsv')
        with self.assertRaises(ParsingError):
            list(parser.parse())

    def test_geostas_parser(self):
        self._create_pdb_file('geostas.pdb', atom_count=6)
        self._create_geostas_clustering_file('clustering_hier_02.json', [
            [1, 2, 3],
            [4, 5, 6],
        ])

        parser = geostas.Parser(self.root_path / 'geostas.pdb', self.root_dir.name)
        result = list(parser.parse())

        self.assertEqual(
            result,
            [('GeoStaS Hierarchical', 2, '1-3,4-6')],
        )

        self._create_geostas_clustering_file('clustering_kmeans_03.json', [
            [1, 2],
            [5],
            [3, 4, 6],
        ])
        result = list(parser.parse())

        self.assertEqual(
            result,
            [
                ('GeoStaS K-means', 3, '1-2,5-5,3-4_6-6'),
                ('GeoStaS Hierarchical', 2, '1-3,4-6'),
            ],
        )

    def test_geostas_chopping(self):
        parser = geostas.Parser('unused.pdb', 'unused/')

        residue_groups = [[1, 2, 3], [10, 11, 20, 21]]
        chops = parser._generate_chopping(residue_groups)

        self.assertEqual(chops, '1-3,10-11_20-21')

        residue_groups = [[1, 2, 3], [5], [10, 11, 20, 21]]
        chops = parser._generate_chopping(residue_groups)

        self.assertEqual(chops, '1-3,5-5,10-11_20-21')

    def _create_pdb_file(self, name, atom_count):
        coordinates = np.zeros((atom_count, 3))

        trajectory = Trajectory.from_ca_frames([coordinates], topology_attr={
            'resids': np.arange(atom_count) + 1,
        })

        trajectory.write_static(self.root_path / name)

    def _create_chainsaw_clustering_file(self, name, data):
        with open(self.root_path / name, 'w') as f:
            writer = csv.writer(f, dialect='unix', delimiter='\t')
            writer.writerow(['chain_id', 'sequence_md5', 'nres', 'ndom', 'chopping', 'confidence', 'time_sec'])
            writer.writerow(data)

    def _create_merizo_clustering_file(self, name, data):
        with open(self.root_path / name, 'w') as f:
            writer = csv.writer(f, dialect='unix', delimiter='\t')
            writer.writerow(data)

    def _create_geostas_clustering_file(self, name, data):
        json_data = json.dumps(data)
        with open(self.root_path / name, 'w') as f:
            f.write(json_data)


if __name__ == '__main__':
    unittest.main()
