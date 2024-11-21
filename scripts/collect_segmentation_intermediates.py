"""
Input:
- <chainsaw/protein.tsv>: Segmentation by chainsaw
- <merizo/protein.tsv>:   Segmentation by merizo
- ...

Output:
- <segmentation.tsv>: A single tsv with all the different segmentations
"""

import sys
import os
import csv
import json
import re
from pathlib import Path

if len(sys.argv) < 2:
    print(f"USAGE: python {os.path.basename(__file__)} [inputs] <output>")
    sys.exit(1)

def main():
    input_files = sys.argv[1:-1]
    output_file = sys.argv[-1]

    columns = ['index', 'method', 'domain_count', 'chopping']
    segmentations = []
    index = 1

    for path in input_files:
        if '/chainsaw/' in path:
            rows = read_csv_rows(path, delimiter='\t')
            data = rows[0]
            segmentations.append((index, "chainsaw", data['ndom'], data['chopping']))
            index += 1
        elif '/merizo/' in path:
            pass
            # Ignore for now, fails on these PDBs
        elif '/bio3d_geostas/' in path:
            # Path is a directory with clustering_kmeans_NN.json files:
            for file in sorted(Path(path).glob('clustering_kmeans_*.json')):
                k           = int(re.findall(r'\d\d', Path(file).stem)[0])
                atom_groups = json.loads(Path(file).read_text())
                chopping    = generate_geostas_chopping(atom_groups)
                description = f"GeoStaS K-means, K={k}"

                segmentations.append((index, description, len(atom_groups), chopping))
                index += 1

            for file in sorted(Path(path).glob('clustering_hier_*.json')):
                k           = int(re.findall(r'\d\d', Path(file).stem)[0])
                atom_groups = json.loads(Path(file).read_text())
                chopping    = generate_geostas_chopping(atom_groups)
                description = f"GeoStaS Hierarchical, K={k}"

                segmentations.append((index, description, len(atom_groups), chopping))
                index += 1

        else:
            print(f"ERROR: Don't know how to process path: {path}")
            sys.exit(1)

    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t', dialect='unix', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(columns)
        for segmentation in segmentations:
            writer.writerow(segmentation)


def read_csv_rows(path, **kwargs):
    with open(path) as f:
        reader = csv.DictReader(f, **kwargs)
        return [row for row in reader]


def generate_geostas_chopping(atom_groups):
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

if __name__ == "__main__":
    main()
