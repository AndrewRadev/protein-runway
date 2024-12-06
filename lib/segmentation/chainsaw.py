import csv

from . import SegmentationParser


class Parser(SegmentationParser):
    def parse(self):
        rows = self._read_csv_rows(self.path, delimiter='\t')
        data = rows[0]

        yield ("Chainsaw", data['ndom'], data['chopping'])

    def _read_csv_rows(self, path, **kwargs):
        with open(self.path) as f:
            reader = csv.DictReader(f, **kwargs)
            return [row for row in reader]
