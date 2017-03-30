""" tests for mpinterfaces.transformations """

from __future__ import unicode_literals

import unittest
from mpinterfaces.mpint_parser import mpint_parse_arguments
from io import StringIO
import sys

__author__ = "Seve G. Monahan"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 28, 2017"


class TestGetRList(unittest.TestCase):

    def setUp(self):
        self.old_stdout = sys.stdout
        sys.stdout = StringIO()

    def tearDown(self):
        sys.stdout = self.old_stdout

    def test_basic(self):
        with self.assertRaises(SystemExit) as cm:
            mpint_parse_arguments(['-h'])

        self.assertEqual(cm.exception.code, 0);

    def test_input(self):
        args = mpint_parse_arguments(['--input', 'test_input.yaml'])
        self.assertEqual(args.input, 'test_input.yaml')


if __name__ == "__main__":
    unittest.main()
