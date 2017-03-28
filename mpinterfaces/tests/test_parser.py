""" tests for mpinterfaces.transformations """

from __future__ import unicode_literals

import unittest
from mpinterfaces.mpint_parser import create_parser

__author__ = "Seve G. Monahan"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 28, 2017"


class TestGetRList(unittest.TestCase):

    def setUp(self):
        self.parser = create_parser()

    def test_basic(self):
        string = self.parser.parse_args(['-h'])
        self.assertTrue(string != "")


if __name__ == "__main__":
    unittest.main()
