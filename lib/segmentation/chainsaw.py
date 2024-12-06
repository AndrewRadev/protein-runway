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

        yield ("Chainsaw", data['ndom'], data['chopping'])


def _read_csv_rows(path, **kwargs):
    with open(path) as f:
        reader = csv.DictReader(f, **kwargs)
        return [row for row in reader]
