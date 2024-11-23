######
# Generates a pymol script to visualize the segmentation from the standardized chopping information
# The residues which are not part of any domain are colored in gray80
# The script can be run from PyMol using the command: run <protein.pml> or by going to File > Run Script

# Input: <protein.tsv> (tsv file from collect_segmentation_intermediates.py)
#        <protein.pdb>
#        <method> (chainsaw, k_geostas, h_geostas)
#        <k> (1-10)
#        <protein.pml> (optional)
# Output: <protein.pml>
######


import sys
import os
from pathlib import Path
import csv

if len(sys.argv) < 3:
    print(f"USAGE: python {os.path.basename(__file__)} <protein.tsv> <protein.pdb> <method> <k> [protein.pml]")
    sys.exit(1)

tsv_file = sys.argv[1]
pdb_file = sys.argv[2]
method = sys.argv[3]
k = sys.argv[4]
output_file = sys.argv[5] if len(sys.argv) > 5 else (f'{Path(pdb_file).stem}.pml')


def parse_chopping(chopping, method, k):
    """
    Columns: 'index', 'method', 'domain_count', 'chopping'
    Example chopping format: '1-100,101-200,201-300_301-400,401-500'
    Method choices: 'chainsaw', 'k_geostas', 'h_geostas'
    k choices: 2-10
    """
    methods = {'chainsaw': 'chainsaw',
               'k_geostas': 'GeoStaS K-means',
               'h_geostas': 'GeoStaS Hierarchical'}
    segments = {}
    with open(tsv_file, 'r') as f:
        # Skip the header
        next(f)
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if methods[method] in row[1] and row[2] == k:
                chopping = row[3]
                break
        if chopping == 'NA':
                print(f"Chopping not available for {method} with k={k}")
                exit(1)
        for i, domain in enumerate(chopping.split(',')):
            r = []
            subdomains = domain.split('_')
            for subdomain in subdomains:
                start, end = subdomain.split('-')
                r.append((int(start) - 1, int(end)))
            segments[i] = r
    return segments


def generate_pymol_script(pdb_file, segmentation, output_file):
    with open(output_file, 'w') as f:
        f.write(f"load {pdb_file}\n")
        f.write(f"color gray80, {Path(pdb_file).stem}\n")
        for i, domain in segmentation.items():
            for j, region in enumerate(domain):
                f.write(f"create domain_{i}_{j}, res {region[0] + 1}-{region[1]}\n")
                f.write(f"color {i + 1}, domain_{i}_{j}\n")
                f.write(f"show cartoon, domain_{i}_{j}\n")
        f.write("hide everything\n")
        f.write("bg_color white\n")
        f.write("zoom all\n")
        f.write("show cartoon, all\n")
        f.write(f"disable {Path(pdb_file).stem}\n")
        f.write(f"enable {Path(pdb_file).stem}\n")


print(f"Generating PyMOL script for {pdb_file}")
print(f"Reading segmentation from {tsv_file}")
segmentation = parse_chopping(tsv_file, method, k)
print(f"Segmentation: {segmentation}")
print(f"Writing PyMOL script to {output_file}")
generate_pymol_script(pdb_file, segmentation, output_file)
exit(0)
