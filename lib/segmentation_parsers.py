from abc import abstractmethod
import csv
import json
import re
from pathlib import Path


def write_segmentations(seg_objects, output_file):
    segmentations = []
    columns = ['index', 'method', 'domain_count', 'chopping']
    index = 1

    for seg_object in seg_objects:
        for method, domain_count, chopping in seg_object.parse():
            segmentations.append((index, method, domain_count, chopping))
            index += 1

    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t', dialect='unix', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(columns)

        for segmentation in segmentations:
            writer.writerow(segmentation)


class SegmentationParser:
    def __init__(self, path):
        self.path = path

    @abstractmethod
    def parse():
        raise NotImplementedError


class ChainsawParser(SegmentationParser):
    def __init__(self, path):
        super().__init__(path)

    def read_csv_rows(self, path, **kwargs):
        with open(self.path) as f:
            reader = csv.DictReader(f, **kwargs)
            return [row for row in reader]

    def parse(self):
        rows = self.read_csv_rows(self.path, delimiter='\t')
        data = rows[0]
        segmentation = ("Chainsaw", data['ndom'], data['chopping'])

        return [segmentation]


class MerizoParser(SegmentationParser):
    def __init__(self, path):
        super().__init__(path)


class GeostasParser(SegmentationParser):
    def __init__(self, path):
        super().__init__(path)

    def generate_geostas_chopping(self, atom_groups):
        """
        Input: [[1, 2, 3], [10, 11, 20, 21], ...]
        Output: 1-3,10-11_20,21,...
        """
        # Structure: { group: [(start, end), (start, end), ...] }
        groupings = {}

        for i, residues in enumerate(atom_groups):
            group = i + 1

            for new_residue in residues:
                if group not in groupings:
                    groupings[group] = [(new_residue, new_residue)]
                    continue

                last_pair = groupings[group][-1]
                start_residue, end_residue = last_pair

                if new_residue - end_residue == 1:
                    # this is the next one, just extend the range:
                    groupings[group][-1] = (start_residue, new_residue)
                else:
                    groupings[group].append((new_residue, new_residue))

        group_descriptions = []
        for i in range(1, len(groupings.keys()) + 1):
            group_descriptions.append('_'.join([f"{start}-{end}" for start, end in groupings[i]]))

        return ','.join(group_descriptions)

    def parse(self):
        segmentations = []

        for file in sorted(Path(self.path).glob('clustering_kmeans_*.json')):
            k           = int(re.findall(r'\d\d', Path(file).stem)[0])
            atom_groups = json.loads(Path(file).read_text())
            chopping    = self.generate_geostas_chopping(atom_groups)
            method      = f"GeoStaS K-means, K={k}"

            segmentations.append((method, len(atom_groups), chopping))

        for file in sorted(Path(self.path).glob('clustering_hier_*.json')):
            k           = int(re.findall(r'\d\d', Path(file).stem)[0])
            atom_groups = json.loads(Path(file).read_text())
            chopping    = self.generate_geostas_chopping(atom_groups)
            method      = f"GeoStaS Hierarchical, K={k}"

            segmentations.append((method, len(atom_groups), chopping))

        return segmentations
