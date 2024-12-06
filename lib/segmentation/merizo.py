import csv
from collections.abc import Iterator
from typing import Tuple

from . import SegmentationParser


class Parser(SegmentationParser):
    def parse(self) -> Iterator[Tuple[str, int, str]]:
        rows = self._read_csv_rows(self.path, delimiter='\t')
        data = rows[0]

        domain_count = data[4]
        chopping     = data[7]

        yield ("Merizo", domain_count, chopping)

    def _read_csv_rows(self, path, **kwargs):
        with open(self.path) as f:
            reader = csv.reader(f, **kwargs)
            return [row for row in reader]
