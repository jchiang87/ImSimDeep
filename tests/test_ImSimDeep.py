"""
Example unit tests for ImSimDeep package
"""
from __future__ import absolute_import, print_function
import os
from collections import namedtuple
import unittest
import desc.imsimdeep

class ObservedSEDsTestCase(unittest.TestCase):
    "TestCase class for ObservedSEDs."
    def setUp(self):
        with open(os.path.join(os.environ['IMSIMDEEP_DIR'], 'tests',
                               'tiny_object_list.txt')) as input_:
            self.object_lines = [x.strip() for x in input_]

    def tearDown(self):
        pass

    def test_get_SED(self):
        "Test the get_SED function."
        observed_seds = desc.imsimdeep.ObservedSEDs()

        sed_names = ('starSED/phoSimMLT/lte037-5.5-1.0a+0.4.BT-Settl.spec.gz',
                     'galaxySED/Inst.32E09.02Z.spec.gz')
        extinction = namedtuple('extinction', 'iA_v iR_v gA_v gR_v'.split())
        ext_pars = (extinction(0, 0, 0.0231692746, 3.1),
                    extinction(0.0273868213, 3.1, 0.5, 3.0999999))
        magnorms = (20.7706993, 17.5603905)
        uband_mags = (23.0524608965, 20.030077675)
        for i, object_line in enumerate(self.object_lines):
            my_sed = observed_seds.get_SED(object_line)
            self.assertTrue(my_sed.name.endswith(sed_names[i]))
            self.assertAlmostEqual(my_sed.magnorm, magnorms[i])
            self.assertAlmostEqual(my_sed.iA_v, ext_pars[i].iA_v)
            self.assertAlmostEqual(my_sed.iR_v, ext_pars[i].iR_v)
            self.assertAlmostEqual(my_sed.gA_v, ext_pars[i].gA_v)
            self.assertAlmostEqual(my_sed.gR_v, ext_pars[i].gR_v)
            self.assertAlmostEqual(my_sed.calcMag('u'), uband_mags[i])

if __name__ == '__main__':
    unittest.main()
