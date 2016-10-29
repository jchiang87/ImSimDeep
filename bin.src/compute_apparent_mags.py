#!/usr/bin/env python
"""
Compute apparent magnitudes for each object in a phosim instance catalog.
"""
from __future__ import absolute_import, print_function
import argparse
import numpy as np
import desc.imsimdeep

parser = argparse.ArgumentParser()
parser.add_argument('instance_catalog', type=str,
                    help='The phosim instance catalog')
parser.add_argument('outfile', type=str, help='The output filename')

args = parser.parse_args()

bands = 'ugrizy'
band = 'r'  # default value

observed_seds = desc.imsimdeep.ObservedSEDs()

output = open(args.outfile, 'w')
with open(args.instance_catalog) as inst_cat:
    for line in inst_cat:
        object_line = line.strip()
        tokens = object_line.split()
        if tokens[0] == 'filter':
            band = bands[int(tokens[1])]
        if tokens[0] != 'object':
            continue
        object_id, ra, dec = tokens[1:4]
        sed = observed_seds.get_SED(object_line)
        output.write('%s  %s  %s  %.8f\n' %
                     (object_id, ra, dec, sed.calcMag(band)))
output.close()
