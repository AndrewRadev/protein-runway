import json
from pathlib import Path
from collections.abc import Iterator
from typing import Tuple

import MDAnalysis as mda

from . import SegmentationParser


class Parser(SegmentationParser):
    def __init__(self, pdb_path, clustering_directory_path):
        super().__init__(pdb_path, clustering_directory_path)

    def parse(self) -> Iterator[Tuple[str, int, str]]:
        pdb_path, clustering_directory_path = self.paths

        mda_universe = mda.Universe(pdb_path)
        mda_atoms    = mda_universe.select_atoms('name = CA')

        for file in sorted(Path(clustering_directory_path).glob('clustering_kmeans_*.json')):
            atom_groups = json.loads(Path(file).read_text())
            chopping    = self._generate_chopping(mda_atoms, atom_groups)
            method      = "GeoStaS K-means"

            yield (method, len(atom_groups), chopping)

        for file in sorted(Path(clustering_directory_path).glob('clustering_hier_*.json')):
            atom_groups = json.loads(Path(file).read_text())
            chopping    = self._generate_chopping(mda_atoms, atom_groups)
            method      = "GeoStaS Hierarchical"

            yield (method, len(atom_groups), chopping)

    def _generate_chopping(self, mda_atoms, atom_groups):
        """
        Input: [[1, 2, 3], [10, 11, 20, 21], ...]
        Output: 1-3,10-11_20,21,...

        Note that the input is (alpha carbon) atom indices (1-indexed) while
        the output is residues. We need to translate the sequential atom index
        into the residue it corresponds to by using an MDAnalysis universe.
        """
        # Structure: { group: [(start, end), (start, end), ...] }
        groupings = {}

        for i, atom_indices in enumerate(atom_groups):
            group = i + 1

            for atom_index in atom_indices:
                new_residue = mda_atoms[atom_index - 1].resnum

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
