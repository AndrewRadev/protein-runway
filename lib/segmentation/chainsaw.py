import csv

from . import SegmentationParser


class Parser(SegmentationParser):
    def __init__(self, path):
        super().__init__(path)

    def read_csv_rows(self, path, **kwargs):
        with open(self.path) as f:
            reader = csv.DictReader(f, **kwargs)
            return [row for row in reader]

    def parse(self):
        rows = self.read_csv_rows(self.path, delimiter='\t')
        data = rows[0]
        segmentation = ("Chainsaw", data['ndom'], data['chopping'])

        return [segmentation]
