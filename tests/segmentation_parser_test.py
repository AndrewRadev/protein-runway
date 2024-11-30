import unittest
import tempfile
import csv
import json
from pathlib import Path

from lib.segmentation_parsers import GeostasParser, ChainsawParser

class TestSegmentationParser(unittest.TestCase):
    def setUp(self):
        self.root_dir = tempfile.TemporaryDirectory()
        self.root_path = Path(self.root_dir.name)

    def tearDown(self):
        self.root_dir.cleanup()

    def test_geostas_chopping(self):
        geostas = GeostasParser('unused/')

        atom_groups = [[1, 2, 3], [10, 11, 20, 21]]
        chops = geostas.generate_geostas_chopping(atom_groups)

        self.assertEqual(chops, '1-3,10-11_20-21')

        atom_groups = [[1, 2, 3], [5], [10, 11, 20, 21]]
        chops = geostas.generate_geostas_chopping(atom_groups)

        self.assertEqual(chops, '1-3,5-5,10-11_20-21')

    def test_chainsaw_parser(self):
        self._create_chainsaw_clustering_file(
            'chainsaw.tsv',
            ['<unused>', '<unused>', 6, 2, '1-3,10-12', '<unused>', '<unused>']
        )

        chainsaw = ChainsawParser(self.root_path / 'chainsaw.tsv')

        result = chainsaw.parse()
        exp_result = [(
            'Chainsaw',
            '2',
            '1-3,10-12',
        )]

        self.assertEqual(result, exp_result)

    def test_geostas_parser(self):
        self._create_geostas_clustering_file('clustering_hier_02.json', [
            [1, 2, 3],
            [10, 11, 12, 50, 60],
        ])

        geostas = GeostasParser(self.root_dir.name)
        result = geostas.parse()
        self.assertEqual(
            result,
            [('GeoStaS Hierarchical, K=2', 2, '1-3,10-12_50-50_60-60')],
        )

        self._create_geostas_clustering_file('clustering_kmeans_03.json', [
            [1, 2, 3],
            [10, 11, 12],
            [50, 60]
        ])
        geostas = GeostasParser(self.root_dir.name)
        result = geostas.parse()

        self.assertEqual(
            result,
            [
                ('GeoStaS K-means, K=3', 3, '1-3,10-12,50-50_60-60'),
                ('GeoStaS Hierarchical, K=2', 2, '1-3,10-12_50-50_60-60'),
            ],
        )

    def _create_chainsaw_clustering_file(self, name, data):
        with open(self.root_path / name, 'w') as f:
            writer = csv.writer(f, dialect='unix', delimiter='\t')
            writer.writerow(['chain_id', 'sequence_md5', 'nres', 'ndom', 'chopping', 'confidence', "time_sec"])
            writer.writerow(data)

    def _create_geostas_clustering_file(self, name, data):
        json_data = json.dumps(data)
        with open(self.root_path / name, 'w') as f:
            f.write(json_data)


if __name__ == '__main__':
    unittest.main()
