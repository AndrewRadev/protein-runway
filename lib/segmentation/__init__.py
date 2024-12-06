from abc import abstractmethod
import csv


def write_segmentations(seg_objects, output_file):
    segmentations = []
    columns = ['index', 'method', 'domain_count', 'chopping']
    index = 1

    for seg_object in seg_objects:
        for method, domain_count, chopping in seg_object.parse():
            segmentations.append((index, method, domain_count, chopping))
            index += 1

    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t', dialect='unix', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(columns)

        for segmentation in segmentations:
            writer.writerow(segmentation)


class SegmentationParser:
    def __init__(self, path):
        self.path = path

    @abstractmethod
    def parse():
        raise NotImplementedError
