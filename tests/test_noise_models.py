#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
import unittest

from obstools_ipgp.noise_models import PetersonNoiseModel, PressureNoiseModel, two_pole_HP


class TestMethods(unittest.TestCase):
    """
    """
#     def setUp(self):
#         self.path = Path(inspect.getfile(
#             inspect.currentframe())).resolve().parent
#         self.test_path = self.path / "data"

    def test_PetersonNoiseModel(self):
        """
        """
        # Test specifically provided values
        # Test interpolation
        low, high = PetersonNoiseModel([0.01, 0.10, 6.0, 70, 600])
        self.assertListEqual(list(low), [-140, -168, -149, -187.5, -184.4])
        low, high = PetersonNoiseModel([0.01, 0.10, 3.8, 15.4, 354.8])
        self.assertListEqual(list(high), [-90, -91.5, -98, -120, -126])
        low, high = PetersonNoiseModel([0.1, 1, 10, 100, 1000])
        self.assertListEqual([int(x) for x in low], [-168, -166, -163, -185, -178])
        self.assertListEqual([int(x) for x in high], [-91, -116, -115, -131, -111])

    def test_PressureNoiseModel(self):
        """
        """
        # Test specifically provided values
        low, high = PressureNoiseModel([0.01, 0.10, 1.0, 10, 100])
        self.assertListEqual(list(low), [64-120, 71-120, 95-120, 114-120, 130-120])
        self.assertListEqual(list(high), [82-120, 95-120, 118-120, 140-120, 170-120])
        # Test interpolation
        low, high = PressureNoiseModel([0.1, 1, 10, 100, 1000])
        self.assertListEqual([int(x) for x in low], [-49, -25, -6, 10, 30])
        self.assertListEqual([int(x) for x in high], [-25, -2, 20, 50, 60])


#     def test_two_pole_HP(self):
#         """
#         """
#         result = two_pole_HP([0.1, 1, 10, 100, 1000], 10, [1, 1, 1, 1, 1])
#         self.assertListEqual(result, [-49, -25, -6, 10, 30])


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
