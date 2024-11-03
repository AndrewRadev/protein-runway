import sys
import os
import csv
from pathlib import Path

import numpy as np
import scipy.cluster.hierarchy as spc

if len(sys.argv) < 3:
    print(f"USAGE: python {os.path.basename(__file__)} <correlation.txt> <output.tsv>")
    sys.exit(1)

correlation_file = sys.argv[1]
output_file      = sys.argv[2]

protein_name = Path(correlation_file).stem

correlation_matrix = np.loadtxt(correlation_file)

# Code adapted from: https://stackoverflow.com/a/76042335
pdist_uncondensed = 1.0 - correlation_matrix
pdist_condensed = np.concatenate([row[i+1:] for i, row in enumerate(pdist_uncondensed)])
linkage = spc.linkage(pdist_condensed, method='complete')
clusters = spc.fcluster(linkage, 0.5 * pdist_condensed.max(), 'distance')

# Structure: { group: [(start, end), (start, end), ...] }
groupings = {}

for i, group in enumerate(clusters):
    new_residue = i + 1

    if group not in groupings:
        groupings[group] = [(new_residue,)]
        continue

    last_pair = groupings[group][-1]

    if len(last_pair) == 1:
        # We only have the first item in the (start, end) pair, add the second one
        start_residue,  = last_pair
        groupings[group][-1] = (start_residue, new_residue)
    elif len(last_pair) == 2:
        start_residue, end_residue = last_pair

        if new_residue - end_residue == 1:
            # this is the next one, just extend the range:
            groupings[group][-1] = (start_residue, new_residue)
        else:
            groupings[group].append((new_residue,))

group_descriptions = []
for i in range(1, max(groupings.keys())):
    group_descriptions.append('_'.join([f"{start}-{end}" for start, end in groupings[i]]))

description = ','.join(group_descriptions)

with open(output_file, 'w') as out:
    writer = csv.writer(out, delimiter='\t', quoting=csv.QUOTE_MINIMAL, dialect='unix')
    writer.writerow(("chain_id", "nres", "ndom", "chopping"))
    writer.writerow((protein_name, len(clusters), max(clusters), description))
