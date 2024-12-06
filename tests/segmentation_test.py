import unittest
import tempfile
import csv
import json
import warnings
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader as MDAMemoryReader
import numpy as np

from lib.segmentation import (
    geostas,
    chainsaw
)

class TestSegmentationParser(unittest.TestCase):
    def setUp(self):
        self.root_dir = tempfile.TemporaryDirectory()
        self.root_path = Path(self.root_dir.name)
        warnings.simplefilter("ignore")

    def tearDown(self):
        self.root_dir.cleanup()
        warnings.resetwarnings()

    def test_geostas_chopping(self):
        parser = geostas.Parser('unused.pdb', 'unused/')

        residue_groups = [[1, 2, 3], [10, 11, 20, 21]]
        chops = parser._generate_chopping(residue_groups)

        self.assertEqual(chops, '1-3,10-11_20-21')

        residue_groups = [[1, 2, 3], [5], [10, 11, 20, 21]]
        chops = parser._generate_chopping(residue_groups)

        self.assertEqual(chops, '1-3,5-5,10-11_20-21')

    def test_chainsaw_parser(self):
        self._create_chainsaw_clustering_file(
            'chainsaw.tsv',
            ['<unused>', '<unused>', 6, 2, '1-3,10-12', '<unused>', '<unused>']
        )

        parser = chainsaw.Parser(self.root_path / 'chainsaw.tsv')
        result = list(parser.parse())

        exp_result = [(
            'Chainsaw',
            '2',
            '1-3,10-12',
        )]

        self.assertEqual(result, exp_result)

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

    def _create_pdb_file(self, name, atom_count):
        coordinates = np.zeros((atom_count, 3))

        u = mda.Universe.empty(
            n_atoms=atom_count,
            n_residues=atom_count,
            n_frames=1,
            atom_resindex=np.arange(atom_count),
        )
        u.load_new(coordinates, format=MDAMemoryReader)
        u.add_TopologyAttr('names', ['CA'] * atom_count)
        u.add_TopologyAttr('resids', np.arange(atom_count) + 1)

        with warnings.catch_warnings(action="ignore"):
            u.atoms.write(self.root_path / name)

    def _create_chainsaw_clustering_file(self, name, data):
        with open(self.root_path / name, 'w') as f:
            writer = csv.writer(f, dialect='unix', delimiter='\t')
            writer.writerow(['chain_id', 'sequence_md5', 'nres', 'ndom', 'chopping', 'confidence', 'time_sec'])
            writer.writerow(data)

    def _create_geostas_clustering_file(self, name, data):
        json_data = json.dumps(data)
        with open(self.root_path / name, 'w') as f:
            f.write(json_data)


if __name__ == '__main__':
    unittest.main()
