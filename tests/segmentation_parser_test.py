import sys
import os
import unittest
from lib.segmentation_parsers import *


class TestSegmentationParser(unittest.TestCase):

    def setUp(self):
        self.chainsaw = ChainsawParser('02_intermediate/chainsaw/1fuu_noPTM.tsv')
        self.geostas = GeostasParser()
    def test_parsers(self):
        pass