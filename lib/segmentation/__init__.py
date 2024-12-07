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
    def __init__(self, *paths: list[Path|str]):
        """
        The intended inputs are filesystem paths where segmentation data can be
        found.

        These could be PDBs or trajectories, but they will most likely be the
        final results of external tools that need to be parsed into a unified
        format.
        """
        self.paths = paths

    @abstractmethod
    def parse(self) -> Iterator[Tuple[str, int, str]]:
        """
        The return value is a collection of items of the form:

            ("<method name>", <number of domains>, "<chopping>")

        They could be collected into a list or yield-ed. The
        `write_segmentations` function only expects that this method returns an
        iterable object.
        """
        raise NotImplementedError


class ParsingError(Exception):
    pass
