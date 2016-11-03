"""
Unit tests for instcat_tools
"""
from __future__ import print_function, absolute_import
import os
from collections import namedtuple
import unittest
import numpy as np
import desc.imsimdeep

class InstcatToolsTestCase(unittest.TestCase):
    "TestCase class for instcat_tools."

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @staticmethod
    def _read_instcat(infile):
        obj = namedtuple('obj', 'objectID ra dec'.split())
        my_objects = []
        with open(infile) as input_:
            for line in input_:
                if not line.startswith('object'):
                    continue
                tokens = line.strip().split()
                my_objects.append(obj(int(tokens[1]), float(tokens[2]),
                                      float(tokens[3])))
        return my_objects

    def test_ang_sep(self):
        "Test ang_sep function."
        ra0, dec0 = 53.0449009, -27.3220807
        for sep in np.logspace(np.log10(0.2/3600.), np.log10(179), 100):
            ra1, dec1 = ra0, dec0 + sep
            if np.abs(dec1) <= 90:
                ang_sep = desc.imsimdeep.ang_sep(ra0, dec0, ra1, dec1)
                self.assertAlmostEqual(ang_sep, sep)

    def test_sky_cone_select(self):
        "Test the sky cone selection code."
        infile = os.path.join(os.environ['IMSIMDEEP_DIR'], 'tests',
                              'tiny_instcat.txt')
        outfile = 'sky_cone_select_output.txt'
        ra = 53.0449009
        dec = -27.3220807
        radius = 0.1
        desc.imsimdeep.sky_cone_select(infile, ra, dec, radius, outfile)
        objs = self._read_instcat(outfile)
        self.assertEqual(len(objs), 1)
        self.assertEqual(objs[0].objectID, 992886536196)
        os.remove(outfile)

if __name__ == '__main__':
    unittest.main()
