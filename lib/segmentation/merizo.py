"""
Merizo predicts domain segmentation based on a deep learning model. Currently,
we're using our own fork to fix an issue with recognizing certain residue codes
like HIP and CYX.

Original source: <https://github.com/psipred/Merizo>
Fork source:     <https://github.com/AndrewRadev/Merizo>
"""

import csv
from collections.abc import Iterator
from typing import Tuple

from . import SegmentationParser, ParsingError


class Parser(SegmentationParser):
    def __init__(self, csv_path):
        super().__init__(csv_path)

    def parse(self) -> Iterator[Tuple[str, int, str]]:
        (csv_path,) = self.paths

        rows = _read_csv_rows(csv_path, delimiter='\t')
        data = rows[0]

        if len(data) < 8:
            raise ParsingError("[Merizo] Segmentation failure")

        domain_count = data[4]
        chopping     = data[7]

        yield ("Merizo", domain_count, chopping)


def _read_csv_rows(path, **kwargs):
    try:
        with open(path) as f:
            reader = csv.reader(f, **kwargs)
            return [row for row in reader]
    except Exception as e:
        raise ParsingError(f"[Merizo] Couldn't parse input path {path}") from e
