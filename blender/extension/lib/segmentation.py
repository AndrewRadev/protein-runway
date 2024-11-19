import csv

def parse_segmentation_file(path):
    """
    Expected columns: 'index', 'method', 'domain_count', 'chopping'
    """
    segmentations = {}

    # TODO (2024-11-17) Handle errors

    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            segmentations[row['method']] = row['chopping']

    return segmentations


def generate_domain_ranges(chopping):
    domain_regions = []

    for domain in chopping.split(','):
        regions = []
        for region in domain.split('_'):
            start, end = region.split('-')

            # The start in the chopping is 1-indexed, so translating to a range means start is -1.
            start_index = int(start) - 1

            # Ranges are end-exclusive, so we do not remove -1 from the end
            end_index = int(end)

            regions.append(range(start_index, end_index))
        domain_regions.append(regions)

    return domain_regions
