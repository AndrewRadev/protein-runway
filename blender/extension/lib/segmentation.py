import csv

def parse_segmentation_file(path):
    """
    Expected columns: 'index', 'method', 'domain_count', 'chopping'
    """
    segmentations = {}

    # TODO (2024-11-17) Handle errors
    # TODO (2024-11-17) Split chopping into lists of ranges

    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            segmentations[row['method']] = row['chopping']

    return segmentations
