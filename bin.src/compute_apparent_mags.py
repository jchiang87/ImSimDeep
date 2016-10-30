#!/usr/bin/env python
"""
Compute apparent magnitudes for each object in a phosim instance catalog.
"""
from __future__ import absolute_import, print_function
import argparse
import desc.imsim
import desc.imsimdeep

parser = argparse.ArgumentParser()
parser.add_argument('instance_catalog', type=str,
                    help='The phosim instance catalog')
parser.add_argument('outfile', type=str, help='The output filename')
parser.add_argument('--numrows', type=int, default=None,
                    help='Number of rows to read from the instance catalog')
args = parser.parse_args()

commands, objs = desc.imsim.parsePhoSimInstanceFile(args.instance_catalog,
                                                    numRows=args.numrows)

band = commands['bandpass']
sed_names = [x[0] for x in objs.groupby('sedName').sedName.unique()]

output = open(args.outfile, 'w')
with open(args.instance_catalog) as inst_cat:
    for sed_name in sed_names:
        app_mag = desc.imsimdeep.ApparentMagnitude(sed_name)
        my_objs = objs.query("sedName=='%s'" % sed_name)
        for i in range(len(my_objs)):
            row = my_objs.iloc[i]
            output.write('%i  %.8f  %.8f  %.8f\n'
                         % (row.objectID, row.ra, row.dec, app_mag(row, band)))
        output.flush()
output.close()
