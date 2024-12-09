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
filename = Path(pdb_file).stem
output_file = sys.argv[5] if len(sys.argv) > 5 else (f'{filename}.pml')


def parse_chopping(chopping, method, k):
    """
    Columns: 'index', 'method', 'domain_count', 'chopping'
    Example chopping format: '1-100,101-200,201-300_301-400,401-500'
    Method choices: 'chainsaw', 'k_geostas', 'h_geostas'
    k choices: 2-10
    """
    methods = {'chainsaw': 'Chainsaw',
               'k_geostas': 'GeoStaS K-means',
               'h_geostas': 'GeoStaS Hierarchical'}
    segments = {}
    with open(tsv_file, 'r') as f:
        # Skip the header
        reader = csv.DictReader(f, delimiter='\t')
        chopping = 'NA'
        for row in reader:
            if methods[method] in row['method'] and row['domain_count'] == k:
                chopping = row['chopping']
                break
        if chopping == 'NA':
            print(f"Chopping not available for {method} with k={k}")
            return segments
        for i, domain in enumerate(chopping.split(',')):
            r = []
            subdomains = domain.split('_')
            for subdomain in subdomains:
                start, end = subdomain.split('-')
                r.append((int(start), int(end)))
            segments[i] = r
    return segments


def generate_pymol_script(pdb_file, segmentation, output_file):
    with open(output_file, 'w') as f:
        f.write(f"load {pdb_file}\n")
        f.write(f"color gray50, {filename}\n")
        for i, domain in segmentation.items():
            # Create selections instead of domains
            for j, region in enumerate(domain):
                f.write(f"select domain_{i:02}_{j:02}, res {region[0]}-{region[1]}\n")
                if i < 5:
                    f.write(f"color {i + 2}, domain_{i:02}_{j:02}\n")
                    # skip color 1 because it's black on black bg
                else:
                    f.write(f"color {i + 3}, domain_{i:02}_{j:02}\n")
                    # skip color 7 because it's the same as 6
                f.write(f"show cartoon, domain_{i:02}_{j:02}\n")
        f.write("hide everything\n")
        f.write("zoom all\n")
        f.write("show cartoon, all\n")
        f.write("deselect\n")


print(f"Generating PyMOL script for {pdb_file}")
print(f"Reading segmentation from {tsv_file}")
segmentation = parse_chopping(tsv_file, method, k)

if not segmentation:
    print("Error: Segmentation information is not available, PyMoL script cannot be generated.")
    exit(1)

print(f"Writing PyMOL script to {output_file}")
generate_pymol_script(pdb_file, segmentation, output_file)
