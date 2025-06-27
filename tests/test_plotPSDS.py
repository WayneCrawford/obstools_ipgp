#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
import unittest

from obstools_ipgp.plotPSDs import _get_shared_nslcs, _fdsn_to_regex


class TestMethods(unittest.TestCase):
    """
    """
#     def setUp(self):
#         self.path = Path(inspect.getfile(
#             inspect.currentframe())).resolve().parent
#         self.test_path = self.path / "data"

    def test___get_shared_nslcs(self):
        """
        """
        self.assertListEqual(sorted(_get_shared_nslcs(['A', 'B', 'C', 'D'],
                                                      ['C', 'D', 'E', 'F'])),
                             ['C', 'D'])

    def test___fdsn_to_regex(self):
        """
        """
        self.assertEqual(_fdsn_to_regex('*'), '^.*$')
        self.assertEqual(_fdsn_to_regex('???'), '^...$')


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
