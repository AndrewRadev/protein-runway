from abc import ABC, abstractmethod
import sys
import os
import csv
import json
import re
from pathlib import Path

class SegmentationParser:
    def __init__(self, path):
        self.path = path
    
    @abstractmethod
    def parse():
        raise NotImplementedError

class ChainsawParser(SegmentationParser):

    def __init__(self, path):
        super().__init__(path)
        self.name = 'chainsaw'

    def read_csv_rows(self, path, **kwargs):
        with open(self.path) as f:
            reader = csv.DictReader(f, **kwargs)
            return [row for row in reader]
    
    def parse(self):
        rows = self.read_csv_rows(self.path, delimiter='\t')
        data = rows[0]
        return data['ndom'], data['chopping']
    
class MerizoParser(SegmentationParser):

    def __init__(self, path):
        super().__init__(path)
        self.name = 'merizo'

class GeostasParser(SegmentationParser):

    def __init__(self, path):
        super().__init__(path)
        self.k = ''
        self.h = ''
        
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
        for file in sorted(Path(self.path).glob('clustering_kmeans_*.json')):
            k           = int(re.findall(r'\d\d', Path(file).stem)[0])
            atom_groups = json.loads(Path(file).read_text())
            chopping    = self.generate_geostas_chopping(atom_groups)
            self.k = f"GeoStaS K-means, K={k}"

            return len(atom_groups), chopping
        
        for file in sorted(Path(self.path).glob('clustering_hier_*.json')):
                k           = int(re.findall(r'\d\d', Path(file).stem)[0])
                atom_groups = json.loads(Path(file).read_text())
                chopping    = self.generate_geostas_chopping(atom_groups)
                self.h = f"GeoStaS Hierarchical, K={k}"

                return len(atom_groups), chopping
                
            