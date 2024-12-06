import json
from pathlib import Path
from collections.abc import Iterator
from typing import Tuple

from . import SegmentationParser


class Parser(SegmentationParser):
    def parse(self) -> Iterator[Tuple[str, int, str]]:
        for file in sorted(Path(self.path).glob('clustering_kmeans_*.json')):
            atom_groups = json.loads(Path(file).read_text())
            chopping    = self._generate_chopping(atom_groups)
            method      = "GeoStaS K-means"

            yield (method, len(atom_groups), chopping)

        for file in sorted(Path(self.path).glob('clustering_hier_*.json')):
            atom_groups = json.loads(Path(file).read_text())
            chopping    = self._generate_chopping(atom_groups)
            method      = "GeoStaS Hierarchical"

            yield (method, len(atom_groups), chopping)

    def _generate_chopping(self, atom_groups):
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
            range_description = [f"{start}-{end}" for start, end in groupings[i]]
            group_descriptions.append('_'.join(range_description))

        return ','.join(group_descriptions)
