import csv

from . import SegmentationParser


class Parser(SegmentationParser):
    def __init__(self, path):
        super().__init__(path)

    def parse(self):
        rows = self._read_csv_rows(self.path, delimiter='\t')
        data = rows[0]

        domain_count = data[4]
        chopping     = data[7]

        segmentation = ("Merizo", domain_count, chopping)

        return [segmentation]

    def _read_csv_rows(self, path, **kwargs):
        with open(self.path) as f:
            reader = csv.reader(f, **kwargs)
            return [row for row in reader]
