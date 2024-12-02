import csv
import re
from collections import defaultdict

def parse_segmentation_file(path):
    """
    Expected columns: 'index', 'method', 'domain_count', 'chopping'
    """
    # A nested dictionary of { <method>: { <domain_count>: <chopping> } }
    segmentations = defaultdict(dict)

    # TODO (2024-11-17) Handle errors

    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            segmentations[row['method']][row['domain_count']] = row['chopping']

    return segmentations


def generate_domain_ranges(chopping):
    domain_regions = []

    for domain in chopping.split(','):
        regions = []
        for region in domain.split('_'):
            if re.match('^\d+$', region) is None:
                start, end = region.split('-')
            else:
                # Only a single domain
                start = region
                end   = region

            start_index = int(start)
            end_index   = int(end)

            regions.append(range(start_index, end_index))
        domain_regions.append(regions)

    return domain_regions
