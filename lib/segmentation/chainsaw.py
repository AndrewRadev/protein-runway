"""
Chainsaw predicts domain segmentation based on a deep learning model. It uses
the ``stride`` tool to determine secondary structure, and produces a TSV file
with a description of the output.

Source: https://github.com/JudeWells/Chainsaw
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

        try:
            data         = rows[0]
            domain_count = data['ndom']
            chopping     = data['chopping']
        except (IndexError, KeyError) as e:
            raise ParsingError("[Chainsaw] Invalid input format") from e

        yield ("Chainsaw", domain_count, chopping)


def _read_csv_rows(path, **kwargs):
    try:
        with open(path) as f:
            reader = csv.DictReader(f, **kwargs)
            return [row for row in reader]
    except Exception as e:
        raise ParsingError(f"[Chainsaw] Couldn't parse input path {path}") from e
