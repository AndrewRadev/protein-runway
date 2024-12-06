import csv
from abc import abstractmethod
from collections.abc import Iterator
from typing import Tuple
from pathlib import Path


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
    def __init__(self, path: Path|str):
        """
        The intended input is a filesystem path where segmentation data can be
        found.

        TODO multiple paths
        """
        self.path = path

    @abstractmethod
    def parse(self) -> Iterator[Tuple[str, int, str]]:
        """
        The return value is ("<method name>", <number of domains>, "<chopping>").

        They could be collected into a list or yield-ed.
        """
        raise NotImplementedError
