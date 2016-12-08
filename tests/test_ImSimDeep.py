"""
Example unit tests for ImSimDeep package
"""
from __future__ import absolute_import, print_function
import os
from collections import namedtuple
import unittest
import pandas as pd
import desc.imsim
import desc.imsimdeep

class ObservedSEDsTestCase(unittest.TestCase):
    "TestCase class for ObservedSEDs."
    def setUp(self):
        infile = os.path.join(os.environ['IMSIMDEEP_DIR'], 'tests',
                              'tiny_instcat.txt')
        self.objects = desc.imsim.parsePhoSimInstanceFile(infile).objects

    def tearDown(self):
        pass

    def test_get_SED(self):
        "Test the get_SED function."
        # Expected u-band magnitudes:
        uband_mags = (23.0524608965, 20.0988206)

        for i in range(len(self.objects)):
            pars = self.objects.iloc[i]
            app_mag = desc.imsimdeep.ApparentMagnitude(pars.sedFilepath)
            self.assertAlmostEqual(app_mag(pars, 'u'), uband_mags[i])

if __name__ == '__main__':
    unittest.main()
