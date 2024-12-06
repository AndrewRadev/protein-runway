import csv
from collections.abc import Iterator
from typing import Tuple

from . import SegmentationParser


class Parser(SegmentationParser):
    def __init__(self, csv_path):
        super().__init__(csv_path)

    def parse(self) -> Iterator[Tuple[str, int, str]]:
        csv_path = self.paths[0]

        rows = _read_csv_rows(csv_path, delimiter='\t')
        data = rows[0]

        domain_count = data[4]
        chopping     = data[7]

        yield ("Merizo", domain_count, chopping)


def _read_csv_rows(path, **kwargs):
    with open(path) as f:
        reader = csv.reader(f, **kwargs)
        return [row for row in reader]
