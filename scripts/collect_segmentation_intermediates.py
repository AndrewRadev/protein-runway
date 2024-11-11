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

if len(sys.argv) < 2:
    print(f"USAGE: python {os.path.basename(__file__)} [inputs] <output>")
    sys.exit(1)

def main():
    output_file = sys.argv[-1]
    input_files = sys.argv[1:-1]

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

if __name__ == "__main__":
    main()
